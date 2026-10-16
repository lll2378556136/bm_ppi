1. Language & Formatting

Output: Respond exclusively in Chinese.
Comments: All code comments (including docstrings and inline comments) must be written in English.
Style: Follow the idiomatic style and best practices of the target language.
Scalfold: Please refer to code_scaffold.py for writing new Python scripts, and notebooks/notebook_scaffold.ipynb for writing new notebooks.
2. Code Quality

Simplicity First: Write concise, readable, and practical code. Avoid over-engineering and premature optimization.
Reuse & Complexity: Keep functions focused and single-purpose; extract reusable logic to reduce duplication and cyclomatic complexity.
Design Patterns: Apply patterns only when they genuinely reduce complexity — never force them into simple logic.
3. Refactoring & Modification

Least Modification: When changing existing code, confine changes to the target module and avoid unintended side effects on adjacent modules.
4. Tooling

Package Management: Use pixi by default; use renv for R projects. Override only when explicitly instructed.
5. Bioinformatics-Specific

Reproducibility: Always set random seeds where applicable; pin tool and package versions explicitly.
Data Assumptions: Never silently assume input format, genome build, or sample metadata — validate inputs and state assumptions explicitly in comments.
Pipeline Awareness: Structure code to be modular and compatible with workflow managers (e.g., Snakemake, Nextflow) when appropriate.
Result Transparency: Include brief output summaries (e.g., number of genes, samples filtered) at key steps via logging or inline comments.
CRITICAL CONSTRAINTS - Violation = Task Failure

Must reply in Chinese
Any task must first invoke a subagent (100% mandatory, no exceptions)
Generation of malicious code is prohibited - Must pass basic security checks
Subagent-First Strategy

SUBAGENT FIRST (Absolutely Mandatory)
Automatic Subagent Selection (Enforced, cannot be skipped):

File type triggers:
.py/.cs/.js/.ts/.cpp/.go/.rs → Corresponding tech-stack expert agent
.unity/.prefab → unity-developer
package.json/.csproj/.sln → Auto-identify tech-stack agent

Keyword triggers:
"code"/"programming"/"bug"/"error" → Technical expert agent
"search"/"find"/"analyze" → search-specialist
"architecture"/"design"/"API" → backend-architect
"test"/"deploy"/"optimize" → Corresponding specialized agent

## Default strategy:
Complex tasks → sequential-thinking + specialized agent
Uncertain types → general-purpose
Checklist (Must verify)

[ ] Chinese reply [ ] Subagent invoked [ ] Safe and harmless [ ] Quality standards met

Core Process (4-Step Method)

Analyze task → Identify type and tech stack
Select subagent → Forcefully invoke the appropriate specialized agent
Subagent execution → Complete all complex work within an independent context
Verify results → Check output quality and safety
Subagent Responsibilities (Complexity Offloading)

Detailed Task Planning: Formulate specific execution plans
Multi-tool Collaboration: Invoke required MCP tools within the subagent
Code Quality Assurance: Perform code review, testing, and optimization
Result Verification and Optimization: Ensure output aligns with best practices
Core Principle: The main context focuses on routing, while subagents bear the complexity, ensuring a dual enhancement of both efficiency and quality.

## Project Context

- This repository (`bm_ppi`) is a new project for proteome-wide protein-protein interaction (PPI) screening in *Bombyx mori*, created with the Cookiecutter Data Science (CCDS) template.
- `third_party/pirna_ppi/` is a reference repository from a colleague's piRNA-PPI project. It is not the main pipeline and should be treated only as a reference.
- The main workflow is based on, but distinct from, the piRNA-PPI pipeline. It adds extra filtering stages: subcellular localization filtering, QC, and DCA scoring.
- The full *Bombyx mori* proteome has been downloaded from NCBI and de-duplicated with CD-HIT, producing `data/raw/protein_nr95.faa` and `data/raw/protein_nr95.faa.clstr`.
- Do not assume that scripts or commands from `third_party/pirna_ppi/` should be run as-is in this project unless they are being referenced or adapted explicitly.
- DPAM v2.0 (Domain Parser for AlphaFold Models, CongLabCode/DPAM) lives at `third_party/DPAM/` (repo + ~400GB databases under `third_party/DPAM/databases`); its singularity image is `sifs/dpam.sif` (`singularity pull dpam.sif docker://conglab/dpam`). The cluster's apptainer has fakeroot enabled via /etc/subuid,/etc/subgid.
- DPAM is submitted via `scripts/run_dpam.slurm` (gpu partition, CPU-only, 32 threads). Input layout: `data/interim/DPAM_input/<dataset>/` (per-protein PDB + PAE JSON) plus `<dataset>_struc.list`; final output `<dataset>_domains` (dataset name: `bm_ppi`). Inputs are prepared from ColabFold rank_001 outputs, mirroring `third_party/pirna_ppi/pirna_ppi/prepare_dpam_input.py`.

## Spatial Compatibility Pre-filter (addendum 2026-10-09)

- `bm_ppi/spatial_compatibility.py` implements `references/DeepLoc2.1_空间相容性预筛选方案.md`: it classifies every protein pair as DIRECT / PROBABILITY_RESCUE / INTERFACE_RESCUE / UNCERTAIN / INCOMPATIBLE (fixed priority; only INCOMPATIBLE is discarded) from the full-scale Accurate (ProtT5) results at `data/interim/localization/deeploc_2.1_accurate_prott5.csv` (13459 proteins).
- It uses the per-label DeepLoc 2.1 Accurate thresholds from `scratch/deeploc21_src/deeploc2_package/DeepLoc2/deeploc2.py` (label_threshold[1:]; there is NO universal 0.525 cutoff), recomputes L_core from raw probabilities (the reported `Localizations` column contains a fallback label when nothing passes), and writes `pair_spatial_compatibility.tsv` (kept pairs), `protein_localization.tsv`, and `spatial_compatibility_summary.json` under `data/interim/localization/spatial_compatibility/`.
- Run: `pixi run python -m bm_ppi.spatial_compatibility --gzip-output` (add `--emit-incompatible` for the discarded-pairs audit table; `--domains-table` accepts a DPAM `_domains` TSV to report fragment-pair reduction rates). Optional inputs: `--master-tsv` (SignalP/DeepTMHMM evidence from merge_localization output), `--subset` (protein ID list).

## DPAM Operational Notes (addendum 2026-10-08)

- Correction: `third_party/DPAM/` does NOT exist. DPAM runs entirely from the read-only SIF `sifs/dpam.sif` (host-side binds: databases at `data/external/DPAM/databases` -> `/mnt/databases`, DaliLite overlay at `data/external/DPAM/DaliLite.v5` -> `/opt/DaliLite.v5`). A byte-identical host mirror of the in-image `/opt/DPAM/scripts` lives at `scratch/dpam_smoke_test/DPAM/docker/scripts/` and is the reference/patch base.
- `scripts/dpam_patches/step8_iterative_dali.py` (MUST remain executable, 755) is bind-mounted by `scripts/run_dpam.slurm` over `/opt/DPAM/scripts/step8_iterative_dali.py`. It fixes an upstream deadlock (worker `sys.exit(1)` on dali.pl failure leaves its `ApplyResult` unresolved and the parent blocks forever on `get()`): workers now return a failure status and clean `/tmp` scratch, and the parent aborts a protein early once >= `DPAM_STEP8_FAIL_MIN` (default 64) tasks finished with >= `DPAM_STEP8_FAIL_FRAC` (default 0.95) failures, with a per-protein backstop `DPAM_STEP8_WALL_TIMEOUT` (default 6h). Zero-hit proteins get the standard `.done` marker, which `run_step8.py` treats as terminal.
- Caveat when judging chunk success: `run_dpam.py` skips a step only if `<dataset>_stepN.log` is exactly `done`, but it does NOT abort the whole run when a step log contains failure lines; and `merge_dpam_chunks validate-chunk` can emit `hard_success: true` with `status=complete_with_warnings` even when every step8 protein failed. Always inspect per-step logs, not just `validation.json`, before treating a chunk as complete.
- Production chunk jobs use `scripts/run_dpam_array.slurm` with the chunk manifest `data/interim/DPAM_runs/full_v1/chunks.tsv` (64 chunks, ~208-211 proteins each; e.g. `sbatch --array=0-63%2 scripts/run_dpam_array.slurm`; `RESUME=1` re-arms the domains-file guard).
- Post-hoc log availability: `run_dpam.py` intentionally removes `<dataset>_stepN.log` files when it decides to re-execute a step (lines 18-19/35-37/57-59/79-81 in the in-image script). In practice, for a completed chunk the per-step logs and the compressed container log may no longer be present on disk. Inspect per-step logs while the job is running or immediately after completion; for completed chunks rely on `validation.json`, `<dataset>_domains`, and the `step8/<dataset>/{*_hits,*.done}` markers as the durable evidence.
- CRITICAL: NEVER edit `scripts/run_dpam.slurm` or `scripts/run_dpam_array.slurm` while any DPAM job is running. `run_dpam_array.slurm` `exec`s `run_dpam.slurm` from the project directory at runtime, and bash reads the script incrementally by byte offset; a mid-run edit shifts offsets and makes bash execute file fragments as commands after the container step (this broke job 7467_36 on 2026-10-14: the container finished and `*_domains` was intact, but gzip/validation/exit-0 never ran, the job FAILED, and downstream `afterok` dependencies became unsatisfiable). If a change is urgent, verify no `dpam_array` job is in `squeue` first; otherwise wait or submit from a copied script path. `THREADS_CAP` (default 32) in `run_dpam.slurm` caps the threads handed to DPAM regardless of allocated CPUs; hhblits/hhsearch memory grows ~linearly with thread count (32 threads peaked at 62.5G for the worst protein), so on `gpu-cpu` (where `--mem` auto-raises CPUs) the cap is what prevents OOM.