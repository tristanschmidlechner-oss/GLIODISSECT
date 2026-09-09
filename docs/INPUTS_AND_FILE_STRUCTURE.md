# Inputs and file structure

This guide describes the supported end-to-end directory workflow. It is the recommended starting point for users who have not worked with HypoxiaBOLD before.

## 1. What HypoxiaBOLD expects

HypoxiaBOLD starts after spatial preprocessing. It does not take raw DICOM data and does not perform motion correction, spatial smoothing, anatomical segmentation, tumour segmentation, or registration estimation.

Before running the package, prepare:

- a realigned BOLD time series;
- BOLD volumes spatially smoothed with a 6-mm FWHM Gaussian kernel;
- grey-matter, white-matter, and CSF probability maps resliced into the BOLD/EPI grid;
- the RespirAct EndTidal recording and scanner-event workbook; and
- for patients, an Oncohabitats segmentation whose NIfTI geometry relates it correctly to the resliced T1 image.

The manuscript implementation and packaged directory runner are specific to the validated double-step hypoxia acquisition:

| Setting | Validated value |
|---|---:|
| Repetition time | 1.8 s |
| Number of BOLD volumes | 160 |
| Spatial smoothing | 6-mm FWHM |
| Temporal robust-local-smoothing window | 12 volumes (21.6 s) |
| First hypoxic step | 60 s |
| Inter-step interval | 40 s |
| Second hypoxic step | 60 s |
| Tissue-probability threshold | GM + WM + CSF >= 0.8 |

The runner currently uses TR = 1.8 s and the 60/40/60-s protocol timing directly. It does not infer these values from the NIfTI header. The code determines the number of scans from the supplied images but was validated with 160 volumes. Other protocols require explicit adaptation and revalidation.

## 2. Required directory tree

Use this structure and preserve the spaces and capitalization shown for `Respiract files` and `Events.xlsx`:

```text
Subject/
|-- BOLD/
|   |-- 32sr...00001.nii
|   |-- 32sr...00002.nii
|   `-- ...                              # validated series: 160 volumes
|-- T1/
|   |-- reslicec1s....nii                # grey-matter probability
|   |-- reslicec2s....nii                # white-matter probability
|   |-- reslicec3s....nii                # CSF probability
|   `-- reslices....nii                  # patient workflow only
|-- Respiract files/
|   |-- Events.xlsx
|   `-- ...EndTidal....txt
`-- segmentationoutput/                  # patient workflow only
    `-- oncohabitats/
        `-- native/
            `-- Segmentation.nii
```

An optional `Subject/bBOLD.mat` is used only when `filteringMode` is explicitly set to `'historical_cache'` for reproduction testing. It is ignored during the default fresh analysis.

## 3. File requirements

### BOLD volumes

- Location and selection pattern: `BOLD/32sr*.nii`.
- Supply the already realigned and 6-mm-smoothed volumes used for analysis.
- The validated layout uses one three-dimensional NIfTI file per time point.
- Files are sorted lexicographically. Use zero-padded volume numbers so lexical order equals acquisition order.
- All volumes must have compatible dimensions and NIfTI geometry.
- The runner converts finite zero-valued BOLD samples to missing values before temporal filtering.

### Tissue-probability maps

Exactly one file must match each pattern:

| Tissue | Required pattern |
|---|---|
| Grey matter | `T1/reslicec1s*.nii` |
| White matter | `T1/reslicec2s*.nii` |
| CSF | `T1/reslicec3s*.nii` |

These must be scaled probability maps, normally ranging from 0 to 1. Each map must match the dimensions and affine matrix of the first BOLD volume; the affine tolerance used by the runner is `1e-5`. No registration is estimated when they do not match.

For healthy controls, the analysis mask is `GM + WM + CSF >= 0.8` with finite values in all three maps.

### Patient T1 and segmentation

The patient workflow additionally requires:

- exactly one `T1/reslices*.nii` image already resliced into the BOLD grid; and
- `segmentationoutput/oncohabitats/native/Segmentation.nii`.

The resliced T1 must match the BOLD dimensions and affine. The segmentation is copied into the new output directory and resampled onto the resliced-T1/BOLD grid using nearest-neighbour interpolation and its existing NIfTI geometry. The package does not estimate a new transformation.

Oncohabitats labels 1, 2, and 3 are combined and retained where the resliced T1 intensity is positive. The patient analysis mask is the union of this lesion mask and the tissue-probability mask. Therefore, tumour voxels can be retained even where tissue probabilities are reduced.

### RespirAct EndTidal file

- Exactly one file must match `Respiract files/*EndTidal*.txt`.
- It must contain numeric records with 23 values per record.
- Field 1 is the timestamp in milliseconds.
- Field 6 is end-tidal O2.
- The recording must start before the first BOLD volume and extend beyond the last BOLD volume.

Duplicate timestamps are removed before interpolation. O2 is interpolated to the BOLD acquisition times using shape-preserving cubic interpolation.

### Events workbook

The filename must be exactly `Respiract files/Events.xlsx`. The numeric worksheet must follow the validated RespirAct export:

- column 1: event time in milliseconds;
- column 2: corresponding scanner volume/trigger index;
- first row: skip/offset event;
- intermediate rows: synchronization events used for the linear time fit; and
- final row: sequence-end event.

At least three rows are required so that an intermediate synchronization range exists. The workbook and EndTidal file must describe the same acquisition.

## 4. Installation

```matlab
packageDir = '/path/to/GLIODISSECT';
spmDir = '/path/to/spm12';

restoredefaultpath
addpath(packageDir)
startup
addpath(spmDir)
run(fullfile(packageDir,'tests','run_all_tests.m'))
```

`startup` adds only the publication-package source folders. Add SPM12 separately. Do not add a large legacy MATLAB-script directory, because duplicate function names can change which implementation MATLAB calls.

## 5. Run one healthy control

```matlab
subjectDir = '/path/to/healthy_control';
outputDir = '/path/to/new/healthy_control_output';

maps = run_subject_directory(subjectDir, outputDir, 'healthy_control');
```

## 6. Run one patient

```matlab
subjectDir = '/path/to/patient';
outputDir = '/path/to/new/patient_output';

maps = run_subject_directory(subjectDir, outputDir, 'patient');
```

The output directory must not already exist and must be outside the input subject directory. Input data may remain read-only. The writer refuses to overwrite existing map files.

## 7. Outputs

Both workflows write `timing.mat`, `maps.mat`, and the following NIfTI maps:

| Output file | Meaning | Unit |
|---|---|---|
| `DTP_O2_map_step1_epiMasked.nii` | Delay to peak, step 1 | volumes |
| `DTP_O2_map_step2.nii` | Delay to peak, step 2 | volumes |
| `DTR_O2_map_step1_epiMasked.nii` | Delay to return, step 1 | volumes |
| `DTR_O2_map_step2_epiMasked.nii` | Delay to return, step 2 | volumes |
| `lag_i10_mask.nii` | Corrected arrival delay | volumes |
| `HypoxiaStep1Map_SS_epiMasked.nii` | Steady-state response, step 1 | % BOLD |
| `HypoxiaStep2Map_SS_epiMasked.nii` | Steady-state response, step 2 | % BOLD |
| `HypoxiaStep1Map_SS_O2norm_epiMasked.nii` | O2-normalised steady state, step 1 | % BOLD/mmHg |
| `HypoxiaStep2Map_SS_O2norm_epiMasked.nii` | O2-normalised steady state, step 2 | % BOLD/mmHg |
| `HypoxiaAvg_step1_epiMasked.nii` | Whole-stimulus average, step 1 | % BOLD |
| `HypoxiaAvg_step2_epiMasked.nii` | Whole-stimulus average, step 2 | % BOLD |
| `HypoxiaAvg_step1_O2norm_epiMasked.nii` | O2-normalised average, step 1 | % BOLD/mmHg |
| `HypoxiaAvg_step2_O2norm_epiMasked.nii` | O2-normalised average, step 2 | % BOLD/mmHg |

Patient runs additionally write the copied/resliced segmentation and `segmentation_mask_evidence.mat` inside the output directory.

## 8. Prepared-array workflow

Users with a different directory layout can bypass file discovery by preparing a MATLAB structure and calling `run_healthy_control` or `run_patient`. In this workflow, `bBOLD` must already be temporally filtered.

| Field | Requirement |
|---|---|
| `bBOLD` | Filtered X-by-Y-by-Z-by-T BOLD array |
| `brainLag` | X-by-Y-by-Z initial lag map in volumes |
| `durationRA` | Five protocol-phase durations in volumes |
| `initialGlobalLag` | Initial global lag in volumes |
| `rowScans` | Scan indices, normally `1:T` |
| `shiftedO2` | O2 vector at the corrected lag, or a callback accepting global lag |
| `greyMatter`, `whiteMatter`, `CSF` | Scaled tissue-probability arrays on the BOLD grid |
| `lesionMask` | Patient workflow only; lesion mask on the BOLD grid |

See `examples/example_prepared_arrays.m`.

## 9. G-adjusted normative analysis

The directory runner creates subject-level EPI maps only. It does not perform T1 reslicing, deformation-field estimation, MNI normalization, construction of a healthy-control atlas, or patient deviation-map reconstruction.

The functions under `src/normative` operate on already prepared and mutually aligned arrays:

- `compute_g_scalar` calculates the subject response scalar in native T1 grey-plus-white matter;
- `fit_g_adjusted_model` fits the voxelwise healthy-control model with subject-mask-aware coverage;
- `apply_g_adjusted_model` calculates expected, difference, z, and thresholded deviation maps; and
- `compute_voxelwise_icc` supports paired reproducibility analysis.

Do not mix EPI-masked steady-state maps with the native-T1 grey-plus-white-matter G-scalar pipeline. Normative inputs must use one documented map, mask, and spatial-normalization chain consistently.

## 10. Common errors

| Error identifier | Likely cause |
|---|---|
| `HypoxiaBOLD:MissingBOLD` | No file matches `BOLD/32sr*.nii` |
| `HypoxiaBOLD:AmbiguousInput` | Zero or multiple files match a required single-file pattern |
| `HypoxiaBOLD:Geometry` | Tissue/T1 dimensions or affine do not match BOLD |
| `HypoxiaBOLD:MissingSegmentation` | Patient Oncohabitats segmentation is absent |
| `HypoxiaBOLD:EndTidalFormat` | EndTidal value count is not divisible by 23 |
| `HypoxiaBOLD:TimeCoverage` | EndTidal recording does not span the BOLD acquisition |
| `HypoxiaBOLD:UnsafeOutput` | Output is the input folder or lies inside it |
| `HypoxiaBOLD:ExistingOutput` | The requested output directory already exists |
| `HypoxiaBOLD:SPMUnavailable` | SPM12 is not on the MATLAB path |

If MATLAB reports ambiguous function resolution, start a clean session with `restoredefaultpath`, run this package's `startup`, and then add only SPM12.
