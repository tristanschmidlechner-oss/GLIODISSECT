# Algorithm specification

## Shared temporal pipeline

1. Apply the cohort-specific analysis mask to the filtered BOLD time series.
2. Estimate initial voxelwise response lag by cross-correlation with the O2 trace.
3. Estimate DTP for the first hypoxic step using the 10% and 90% response crossings and bounded iterative baseline/step refinement.
4. Correct the initial lag by the difference between the expected and fitted 10% crossing.
5. Define the global corrected lag as the fifth percentile of finite corrected voxel delays inside the analysis mask.
6. Estimate DTR for step 1, then DTP and DTR for step 2.
7. Define the voxelwise steady-state window from DTP90 to the earlier of DTR90 and the voxel-adjusted stimulus end.
8. Express steady-state signal change as percent BOLD relative to the step-1 baseline and divide by the achieved PetO2 change for the O2-normalised maps.

## Bounded DTP refinement

The package preserves the audited historical DTP update sequence because all 160 archived maps in 20 healthy controls were reproduced from that implementation. The cohort audit found a median of two refinements and criterion-unsatisfied terminal fractions of 2.172% (step 1) and 2.260% (step 2), predominantly two-state cycles. These terminal estimates are therefore described as bounded outputs, not universally converged estimates.

## DTR refinement

DTR independently updates the post-step baseline estimate until the relative baseline-mean change is strictly below 0.0001 or the 15-update bound is reached. All evaluable DTR estimates in the 20-control audit met the criterion within 11 updates.

## G-adjusted normative model

For each voxel, ordinary least squares fits

`steady_state = beta0 + beta1 * G + error`,

where `G` is the subject median O2-normalised steady-state response in the specified native-T1 grey-plus-white-matter mask. The corrected implementation:

- treats each subject mask as voxel-level validity;
- restricts fitting to a scaled anatomical GM+WM mask;
- requires the configured subject coverage (80% by default);
- reports ordinary OLS R-squared without a variance floor;
- applies a fifth-percentile positive residual-SD floor only when standardising patient deviations.

Patient deviations are `z = (observed - expected) / max(raw residual SD, floor)`. R-squared is descriptive and is not used to compute expected, difference or z maps.
