# 基于 DeepLoc 2.1 的蛋白片段空间相容性预筛选方案

**版本：** v1.0  
**日期：** 2026-10-09  
**用途：** 在蛋白片段配对、结构建模或蛋白互作计算之前，利用亚细胞定位信息删除明显缺乏空间相容性的片段组合，从而减少后续计算量。

## 1. 方案结论

本方案将 DeepLoc 2.1 的定位结果用于“空间相容性预筛选”，而不是直接预测蛋白是否互作。

核心规则为：

> 如果两个蛋白具有相同的高可信亚细胞定位标签，则认为它们具有直接空间相容性。

对于最终标签没有交集的蛋白对，不立即全部删除，而是继续检查：

1. 原始定位分数是否提示某个共同区室可能因阈值截断而漏标；
2. 是否存在膜—胞质或膜—胞外等真实空间界面；
3. 定位预测是否存在明显不确定性。

最终将蛋白对分为：

```text
DIRECT             直接空间相容
PROBABILITY_RESCUE 概率容错相容
INTERFACE_RESCUE   界面相容
UNCERTAIN          定位不确定，暂不删除
INCOMPATIBLE       高可信空间不相容，可删除
```

只有 `INCOMPATIBLE` 蛋白对下的片段组合进入删除集合。

## 2. DeepLoc 2.1 输出的使用方式

DeepLoc 2.1 是多标签预测器，可为同一蛋白预测多个亚细胞定位，并额外预测膜关联类型。其定位类别包括：

```text
Nucleus
Cytoplasm
Extracellular
Mitochondrion
Cell membrane
Endoplasmic reticulum
Chloroplast
Golgi apparatus
Lysosome/Vacuole
Peroxisome
```

膜关联类型包括：

```text
Transmembrane
Peripheral
Lipid anchor
Soluble
```

DeepLoc 2.1 官方页面会提供各定位标签的分数和阈值。因此，不建议只保留最高分的一个定位，也不建议只使用未经记录的人工阈值。

参考：

- [DeepLoc 2.1 官方服务说明](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/)
- [DeepLoc 2.1 方法论文](https://academic.oup.com/nar/article/52/W1/W215/7642068)
- [开放获取全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11223819/)

## 3. 数据结构

### 3.1 蛋白级数据

每条完整蛋白至少保存以下字段：

| 字段                  | 含义                    |
| ------------------- | --------------------- |
| `protein_id`        | 蛋白唯一编号                |
| `sequence`          | 完整蛋白序列                |
| `loc_score[k]`      | 区室 `k` 的 DeepLoc 原始分数 |
| `loc_threshold[k]`  | 区室 `k` 的 DeepLoc 判定阈值 |
| `loc_core`          | 超过对应阈值的定位标签集合         |
| `loc_border`        | 接近阈值但未超过阈值的标签集合       |
| `membrane_type`     | 膜关联类型及其分数             |
| `sorting_signals`   | DeepLoc 预测的排序/定位信号    |
| `prediction_model`  | Fast 或 Slow 模式        |
| `prediction_source` | DeepLoc 2.1           |

### 3.2 片段级数据

默认情况下，片段继承完整蛋白的定位标签：

```text
fragment.loc_core       = parent.loc_core
fragment.loc_border     = parent.loc_border
fragment.membrane_type  = parent.membrane_type
fragment.tag_source     = inherited
fragment.parent_id       = protein_id
```

同时保留片段起止坐标：

```text
fragment_id
parent_id
start
end
fragment_sequence
```

这种继承表示“片段属于该蛋白的原始空间背景”，不表示该片段单独表达后一定仍然保持相同定位。

## 4. 定位标签分层

### 4.1 高可信标签

定义：

```text
L_core(P) = { k | p_P(k) >= τ_k }
```

其中：

- `p_P(k)` 是蛋白 `P` 对区室 `k` 的 DeepLoc 分数；
- `τ_k` 是区室 `k` 对应的标签阈值。

### 4.2 近阈值标签

为处理潜在漏标，定义：

```text
L_border(P) = { k | τ_k - δ_k <= p_P(k) < τ_k }
```

`δ_k` 是近阈值缓冲区，建议按标签分别校准。初步测试可以设置一个小的统一缓冲区，但正式分析应尽量使用验证集确定。

### 4.3 不确定状态

以下情况不宜直接判定为不相容：

- `L_core(P)` 为空；
- 只有一个很弱的定位标签；
- 多个标签分数接近；
- 蛋白具有明显多定位或刺激依赖定位特征；
- 片段很短，缺少完整定位信号。

这些蛋白对进入 `UNCERTAIN` 类别。

## 5. 蛋白对的空间相容性判定

### 5.1 直接空间相容：`DIRECT`

两个蛋白的高可信标签有交集：

```text
DIRECT(A, B) = 1, if L_core(A) ∩ L_core(B) != ∅
```

例如：

```text
A = Cytoplasm + Nucleus
B = Nucleus
```

结果：

```text
DIRECT
```

这是最强的空间相容类别，但仍然不代表两个蛋白一定互作。

### 5.2 概率容错相容：`PROBABILITY_RESCUE`

当最终标签没有交集，但同一定位标签在两个蛋白中的原始分数都接近阈值时，可认为可能存在阈值造成的漏标。

将近阈值程度标准化为：

```text
q_P(k) = clip((p_P(k) - (τ_k - δ_k)) / δ_k, 0, 1)
```

其中：

- `q_P(k)=0` 表示分数位于缓冲区下界或更低；
- `q_P(k)=1` 表示分数达到正式阈值；
- `clip` 将结果限制在 0 到 1 之间。

共同定位容错分数定义为：

```text
S_probability(A, B) = max_k min(q_A(k), q_B(k))
```

若：

```text
L_core(A) ∩ L_core(B) = ∅
```

但：

```text
S_probability(A, B) >= θ_probability
```

则标记为：

```text
PROBABILITY_RESCUE
```

`p_A(k) × p_B(k)` 可以作为排序指标，但不应直接解释为真实的共同定位概率。正式阈值 `δ_k` 和 `θ_probability` 应通过已知互作数据校准。

### 5.3 界面相容：`INTERFACE_RESCUE`

膜—胞质关系不应简单视为普通标签交集，也不应一律视为空间不相容。

当满足以下条件时，保留为界面候选：

```text
一方包含 Cell membrane
另一方包含 Cytoplasm
膜蛋白的 membrane_type 为：
Transmembrane、Peripheral 或 Lipid anchor
```

结果标记为：

```text
INTERFACE_RESCUE
```

该类别的含义是：蛋白可能通过膜蛋白的胞质侧结构域、膜结合面或膜锚定区域接触。由于 DeepLoc 2.1 通常不能直接确定膜蛋白的胞质侧、胞外侧或腔面朝向，因此此类结果应降低优先级，而不应与 `DIRECT` 等同。

同样的原则可用于：

```text
Cell membrane + Extracellular
```

但是否保留应结合信号肽、跨膜区和拓扑预测结果。

### 5.4 高可信空间不相容：`INCOMPATIBLE`

建议只有在以下条件同时成立时删除：

```text
1. L_core(A) ∩ L_core(B) = ∅；
2. S_probability(A, B) < θ_probability；
3. 不满足任何界面相容规则；
4. 两个蛋白均有明确的高可信定位；
5. 两个蛋白均不属于明显的多定位或不确定状态。
```

形式化表示为：

```text
INCOMPATIBLE(A, B) =
    no direct overlap
    and no probability rescue
    and no interface rescue
    and high-confidence predictions
```

## 6. 推荐判定优先级

判定顺序建议固定为：

```text
第一步：DIRECT
第二步：PROBABILITY_RESCUE
第三步：INTERFACE_RESCUE
第四步：UNCERTAIN
第五步：INCOMPATIBLE
```

伪代码如下：

```python
def classify_space(A, B):
    # 1. 高可信标签交集
    if set(A.loc_core) & set(B.loc_core):
        return "DIRECT"

    # 2. 原始分数的近阈值容错
    probability_score = common_probability_score(A, B)
    if probability_score >= THETA_PROBABILITY:
        return "PROBABILITY_RESCUE"

    # 3. 膜—胞质或膜—胞外界面
    if membrane_interface(A, B):
        return "INTERFACE_RESCUE"

    # 4. 任一方定位不明确或明显多定位
    if is_uncertain(A) or is_uncertain(B):
        return "UNCERTAIN"

    # 5. 只有明确不相容时才删除
    return "INCOMPATIBLE"
```

## 7. 片段配对流程

空间筛选应在蛋白对层面先进行，再决定是否展开到片段对层面。

```python
for protein_A, protein_B in all_protein_pairs:

    relation = classify_space(protein_A, protein_B)

    if relation == "INCOMPATIBLE":
        continue

    for fragment_A in protein_A.fragments:
        for fragment_B in protein_B.fragments:
            save_candidate(
                fragment_A,
                fragment_B,
                space_class=relation
            )
```

这样，如果两个完整蛋白不相容，就可以一次性删除它们之间的全部片段组合，而不必逐一运行片段级计算。

## 8. 计算量评估

设蛋白 `i` 有 `m_i` 个片段，蛋白 `j` 有 `n_j` 个片段。

未筛选时的片段配对数量：

```text
N0 = Σ(m_i × n_j)
```

筛选后的配对数量：

```text
N1 = Σ(m_i × n_j)，仅统计非 INCOMPATIBLE 蛋白对
```

计算量减少率：

```text
Reduction = 1 - N1 / N0
```

同时记录：

```text
pair_reduction_rate
fragment_pair_reduction_rate
known_interaction_retention
false_exclusion_rate
downstream_runtime
```

## 9. 参数校准与验证

### 9.1 建议的验证数据

优先使用：

- 已知的家蚕蛋白互作对；
- 同一细胞系或相近组织来源的互作数据；
- 有实验亚细胞定位证据的蛋白对；
- 随机匹配或经过表达背景匹配的非互作对作为负集。

### 9.2 需要比较的版本

至少比较三种过滤器：

```text
Baseline：不使用空间筛选
Strict：只使用标签交集
Full：标签交集 + 概率容错 + 界面容错
```

### 9.3 参数选择原则

不要单纯追求最大计算量减少率。建议选择：

```text
在已知互作保留率达到项目要求的前提下，
最大化片段配对减少率。
```

需要重点检查以下亚组：

```text
细胞膜—胞质蛋白对
多定位蛋白对
核—质穿梭蛋白对
分泌/胞外蛋白对
含跨膜片段的蛋白对
```

## 10. 片段级例外处理

完整蛋白标签继承是默认规则，但以下片段应单独标记：

### 10.1 含定位信号的片段

包括：

- 信号肽；
- 线粒体靶向序列；
- 核定位信号；
- 内质网定位或滞留信号；
- 跨膜螺旋；
- 脂质修饰或膜锚定区域。

### 10.2 极短片段

极短片段可能不包含完整定位信号，也不适合单独进行可靠的 DeepLoc 预测。此时保留父蛋白标签，但降低标签置信度，并将来源写为：

```text
tag_source = inherited_parent_context
```

### 10.3 片段单独表达的情况

如果后续实验或计算模型代表的是“片段单独表达”，不能无条件继承父蛋白全部定位标签。建议采用：

```text
父蛋白标签：保留为背景先验
片段预测：作为局部定位证据
两者冲突：标记为 fragment_localization_conflict
```

## 11. 结果文件建议

建议输出三个文件或三个数据表：

### 11.1 `protein_localization.tsv`

保存每个完整蛋白的 DeepLoc 原始结果、阈值、核心标签和近阈值标签。

### 11.2 `fragment_localization.tsv`

保存片段坐标、继承标签、膜关联类型、标签来源和片段级例外标记。

### 11.3 `pair_spatial_compatibility.tsv`

建议字段：

| 字段                         | 含义                                   |
| -------------------------- | ------------------------------------ |
| `protein_A`                | 蛋白A                                  |
| `protein_B`                | 蛋白B                                  |
| `shared_core_labels`       | 共同高可信标签                              |
| `probability_rescue_label` | 概率容错使用的标签                            |
| `interface_type`           | 界面类型                                 |
| `space_class`              | DIRECT/RESCUE/UNCERTAIN/INCOMPATIBLE |
| `space_score`              | 空间相容性排序分数                            |
| `decision`                 | KEEP 或 DISCARD                       |
| `reason`                   | 判定原因                                 |

## 12. 方案边界

本方案只能回答：

> 两个蛋白是否存在合理的空间共存或空间接触可能。

不能单独回答：

- 两个蛋白是否真的结合；
- 是否存在合适的结构界面；
- 是否在同一组织或发育阶段表达；
- 是否具有相同的诱导条件；
- 膜蛋白的具体拓扑朝向；
- 片段单独表达后是否保留完整蛋白定位。

因此，`DIRECT`、`PROBABILITY_RESCUE` 和 `INTERFACE_RESCUE` 都表示“允许进入后续计算”，而不是“确认发生互作”。

## 13. 最终推荐规则

实际运行时采用以下简化版：

```text
1. DeepLoc 2.1 输出多标签和原始分数；
2. 用标签特异性阈值生成 L_core；
3. L_core 有交集：DIRECT，保留；
4. 无交集但共同标签分数接近阈值：PROBABILITY_RESCUE，保留；
5. 无交集但存在有膜关联证据的膜—胞质/膜—胞外关系：INTERFACE_RESCUE，保留；
6. 定位不明确：UNCERTAIN，保留但降低优先级；
7. 只有高可信、无共同定位、无容错关系的蛋白对才判为 INCOMPATIBLE 并删除。
```

该方案保留了“共同定位标签是主要相容依据”的简洁逻辑，同时避免因 DeepLoc 阈值截断、次要定位漏标或膜界面关系而产生过度过滤。
