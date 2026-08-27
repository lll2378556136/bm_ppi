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