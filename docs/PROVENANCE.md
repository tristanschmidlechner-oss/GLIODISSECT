# Source provenance and release scope

The publication package was refactored from the confirmed healthy-control and MATLAB R2025b patient implementations used for the manuscript analyses. It retains one shared DTP/DTR/steady-state numerical core with cohort-specific masking adapters. Its temporal-analysis lineage began within the Respiract Research Group Zurich, now the ASTRAN Lab.

Included source categories are:

- temporal filtering and oxygen-trace preparation;
- initial lag estimation;
- DTP, DTR, corrected-delay, steady-state, and whole-stimulus-average calculations;
- healthy-control and patient mask construction;
- SPM-based NIfTI input/output;
- G-adjusted normative modelling and deviation scoring;
- sensitivity-configuration and voxelwise ICC utilities; and
- synthetic regression and static tests.

Excluded material comprises ROI and heterogeneity analyses, biopsy/DSC workflows, figures, exploratory scripts, duplicate helper versions, participant data, private validation logs, and machine-specific dependency manifests.

The v1.0 numerical core is unchanged from the exactly validated pre-release implementation. The only behavior-affecting code change for public distribution is a platform-neutral output-safety check in `run_subject_directory`: output must be a new directory outside the input subject directory, rather than a machine-specific drive location. This check does not alter filtering, masks, timing, map calculation, or NIfTI encoding. Mojibake in comments and warning messages was also corrected without changing calculations.
