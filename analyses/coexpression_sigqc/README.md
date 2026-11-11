# Co-expression of Axis A members

Runs sigQC on Axis A in eight independent datasets and compares its co-expression with expression-matched random gene sets.

**Used in:** Fig 3b, Table C

| Script | Content |
| --- | --- |
| [`01_build_matrices.R`](code/01_build_matrices.R) | expression matrices |
| [`02_run_sigqc.R`](code/02_run_sigqc.R) | sigQC |
| [`03_summarise.R`](code/03_summarise.R) | summary and random-set comparison |
