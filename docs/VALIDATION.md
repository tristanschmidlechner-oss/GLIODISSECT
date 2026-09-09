# Validation

## Environment

- MATLAB R2025b Update 4
- SPM12
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Regression-validation design

One healthy-control case and one patient case were tested at two levels:

1. **Cached core:** the shared numerical core used the validated historical filtered-BOLD cache.
2. **Fresh R2025b:** filtering was recomputed from prepared spatial inputs using the documented 12-volume robust local smoother.

Each run was compared voxel-by-voxel with an independently produced reference. The patient test also compared the segmentation-derived analysis mask with its reference.

## Results

| Validation level | Healthy control | Patient |
|---|---:|---:|
| Primary maps, cached core | 9/9 exact | 9/9 exact |
| Primary maps, fresh R2025b | 9/9 exact | 9/9 exact |
| Average comparators, cached core | 4/4 exact | 4/4 exact |
| Average comparators, fresh R2025b | 4/4 exact | 4/4 exact |

All comparisons had zero differing values and identical finite/nonfinite patterns, dimensions, affines, datatype, scaling, and zero/NaN encoding. The patient segmentation-derived mask matched exactly. Six unit/static tests passed.

The validated primary maps were DTP steps 1 and 2, DTR steps 1 and 2, corrected delay, steady-state steps 1 and 2, and O2-normalised steady-state steps 1 and 2. The comparator set comprised raw and O2-normalised whole-stimulus averages for both steps.

## Supporting method audits

- **Temporal filter:** among 8-, 12-, and 16-volume robust local-smoothing windows, 12 volumes best reproduced the retained filtered-BOLD reference.
- **Sensitivity:** 20 healthy controls completed the prespecified threshold, temporal-window, and spatial-smoothing configurations.
- **Normative model:** corrected original-map analysis yielded median in-sample OLS R-squared 0.468448 and median leave-one-out Q-squared 0.321003; predictive RMSE improved by 14.56% relative to a mean-only model and improved in 19/20 controls.
- **Iteration behavior:** 160/160 archived arrays across 20 healthy controls reproduced exactly. DTP criterion-unsatisfied fractions had medians of 2.172% for step 1 and 2.260% for step 2; all evaluable DTR estimates met their criterion within 11 refinements.

## Scope

This evidence establishes numerical regression reproduction for the documented implementation and supported layouts. It does not constitute clinical validation, certify unrelated legacy analyses, or guarantee compatibility with arbitrary preprocessing pipelines.
