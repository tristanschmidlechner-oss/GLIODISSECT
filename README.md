# GLIODISSECT: Temporal response decomposition mapping

GLIODISSECT is a MATLAB framework that dissects voxel-wise BOLD time-series responses acquired during a controlled two-step isocapnic hypoxic stimulus into complementary measures of temporal dynamics and steady-state response. Rather than reducing each voxel’s response to a single average value, the framework separately maps response timing, transition dynamics, corrected delay, and steady-state magnitude.

The framework estimates delay-to-peak (DTP), delay-to-return (DTR), corrected delay, steady-state response, oxygen-normalized steady-state response, and whole-stimulus-average comparator maps. Healthy-control and patient workflows use dedicated masking procedures while sharing the same validated numerical core.

GLIODISSECT was developed within the [ASTRAN Lab](https://astranlab.com). This release focuses on the validated temporal decomposition and normative-mapping workflow; downstream ROI, heterogeneity, biopsy, DSC, and unrelated exploratory analyses are outside its scope.

## Requirements

- MATLAB R2025b
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- SPM12 for NIfTI input/output and patient-segmentation reslicing

Curve Fitting Toolbox is not required. Temporal smoothing uses `smoothdata(...,'rloess',12)`.

## Start here

New users should first read [Inputs and file structure](docs/INPUTS_AND_FILE_STRUCTURE.md). It explains what preprocessing must already be complete, the exact healthy-control and patient directory layouts, RespirAct formats, output maps, and common errors.

Choose the entry point that matches your data:

| Starting point | Entry point |
|---|---|
| Preprocessed subject in the supported directory layout | `run_subject_directory` |
| Prepared healthy-control MATLAB arrays | `run_healthy_control` |
| Prepared patient MATLAB arrays | `run_patient` |
| Already aligned maps for normative analysis | Functions under `src/normative` |

## Installation and tests

```matlab
packageDir = '/path/to/GLIODISSECT';
spmDir = '/path/to/spm12';

restoredefaultpath
addpath(packageDir)
startup
addpath(spmDir)
run(fullfile(packageDir,'tests','run_all_tests.m'))
```

## Running prepared arrays

```matlab
input = load('/path/to/prepared_subject.mat');
maps = run_healthy_control(input);   % healthy-control mask
% maps = run_patient(input);         % tissue-plus-lesion mask

save_maps_mat('/path/to/output/subject_maps.mat', maps, ...
    struct('package_version','1.0.0'));
```

Both entry points call `run_hypoxia_bold_core`. See the input contract below and `examples/example_prepared_arrays.m`.

## Running the supported subject-directory layout

```matlab
subjectDir = '/path/to/read_only_subject';
outputDir = '/path/to/new_output_directory';
maps = run_subject_directory(subjectDir, outputDir, 'healthy_control');
% Use 'patient' for the patient segmentation workflow.
```

`outputDir` must not already exist and must be outside `subjectDir`. The runner expects:

- realigned and spatially smoothed `BOLD/32sr*.nii` volumes;
- EPI-space tissue-probability maps under `T1`;
- `Respiract files/Events.xlsx` and an EndTidal text file; and
- for patients, native Oncohabitats `segmentationoutput/oncohabitats/native/Segmentation.nii` and a resliced T1 image.

The exact required filenames and tree are documented in [Inputs and file structure](docs/INPUTS_AND_FILE_STRUCTURE.md). A copy-paste template is provided in `examples/run_subject_directory_example.m`.

The directory runner recomputes temporal filtering and ignores an existing `bBOLD.mat`. It does not estimate motion correction, spatial smoothing, tissue or lesion segmentation, or registration. The patient adapter reslices segmentation labels 1/2/3 by nearest neighbour and retains them where the resliced T1 has positive intensity.

For a filtering-independent reproduction check only, pass `'historical_cache'` as the fourth argument. This explicitly loads the subject-root `bBOLD.mat`.

## Prepared-array input contract

| Field | Description |
|---|---|
| `bBOLD` | Filtered BOLD array, X-by-Y-by-Z-by-T |
| `brainLag` | Initial voxelwise lag in volumes |
| `durationRA` | Five protocol-phase durations in volumes |
| `initialGlobalLag` | Initial global lag in volumes |
| `rowScans` | Scan indices, normally `1:T` |
| `shiftedO2` | O2 trace sampled at corrected lag, or a callback returning it |
| `greyMatter`, `whiteMatter`, `CSF` | Scaled tissue-probability maps in EPI space |
| `lesionMask` | Patient runner only; lesion segmentation in EPI space |

Raw BOLD filtering is available through `filter_bold_timeseries`.

The packaged directory workflow is specific to the manuscript protocol: TR 1.8 s, 160 validated volumes, two 60-s hypoxic steps separated by 40 s, 6-mm spatial smoothing, and a 12-volume temporal smoothing window. These settings are not inferred automatically from arbitrary datasets.

## Outputs

The primary NIfTI outputs comprise DTP for both steps, DTR for both steps, corrected delay, steady-state response for both steps, and O2-normalised steady-state response for both steps. The package also writes raw and O2-normalised whole-stimulus-average comparators for both steps.

Timing maps are expressed in volumes. Response maps are percent BOLD change; O2-normalised maps are percent BOLD change per mmHg.

Output encoding follows the validated historical implementations. Masked DTP step 1, DTR, steady-state, and average-response zeros are written as NaN; corrected-delay zero remains valid. DTP step 2 preserves historical background encoding, so use `maps.analysis_mask` for scientific summaries.

## Methodological notes

- DTP uses bounded refinement. In the 20-control audit, most voxels required two refinements; small terminal subsets reached the update bound without satisfying the post-hoc stability criterion. These values are bounded outputs, not universally converged estimates.
- DTR uses a separate relative baseline-mean stopping criterion of `<0.0001`; all evaluable estimates in the audit met it within 11 refinements.
- The G-adjusted normative model uses ordinary voxelwise OLS. A residual-SD floor is used only to stabilize patient z scores and is not used in R-squared estimation.
- Sensitivity configurations and voxelwise ICC utilities are included for reproducibility analyses.

See [Algorithm specification](docs/ALGORITHM.md), [Validation](docs/VALIDATION.md), and [Provenance](docs/PROVENANCE.md).

## Validation

The frozen numerical implementation was validated in MATLAB R2025b Update 4 using one healthy-control and one patient case at two levels: cached-core reproduction and fresh R2025b filtering. All 36 primary-map comparisons and all 16 whole-stimulus comparator comparisons were exactly equal to independent references, including finite-value patterns, geometry, datatype, and encoding. The patient segmentation-derived mask also matched its independent reference exactly. Six unit/static tests passed.

This is a two-case regression-validation scope; it is not a claim of clinical validation or universal compatibility with arbitrary directory layouts.

## Data and privacy

No participant data, participant identifiers, local filesystem paths, or private validation logs are distributed. Users must obtain and process their own data under appropriate ethics and data-governance approvals.

## Citation and licence

Please cite this software release and the associated methods preprint using the metadata in `CITATION.cff`.

Associated methods preprint: Schmidlechner et al. “Voxel-wise temporal decomposition of hypoxia-targeted BOLD MRI: method development and proof-of-concept application in glioblastoma.” *medRxiv* (2026). [https://doi.org/10.64898/2026.05.27.26354265](https://doi.org/10.64898/2026.05.27.26354265). This preprint has not undergone peer review. After journal publication, cite the peer-reviewed article instead of the preprint.

Software release authors: Tristan Schmidlechner, Vittorio Stumpo, Bas van Niftrik, and Jorn Fierstra; ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich.

Released under the BSD 3-Clause License. See `LICENSE`.
