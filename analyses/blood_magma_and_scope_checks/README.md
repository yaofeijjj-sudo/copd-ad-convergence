# Blood, MAGMA and disease-scope checks

Tests Axis A in peripheral blood with CAMERA and participant-level scores, checks AddNeuroMed duplicates, conditions MAGMA on the Patir signature, restricts the MAGMA background and cross-tabulates the GSE150910 batches.

**Used in:** Fig 6c, Fig 7, Table J, Table K

| Script | Content |
| --- | --- |
| [`01_blood_camera_effects.R`](code/01_blood_camera_effects.R) | blood CAMERA and score differences |
| [`02_addneuromed_followup.R`](code/02_addneuromed_followup.R) | AddNeuroMed duplicate follow-up |
| [`03_magma_background_patir.R`](code/03_magma_background_patir.R) | restricted background and Patir conditioning |
| [`03b_extract_patir.py`](code/03b_extract_patir.py) | Patir core microglia genes |
| [`04_chp_batch_crosstab.R`](code/04_chp_batch_crosstab.R) | GSE150910 diagnosis by batch |
