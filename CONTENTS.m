% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
% GLIODISSECT temporal response decomposition mapping framework, version 1.0.0
%
% Entry points
%   run_subject_directory     - Run a supported preprocessed subject folder.
%   run_healthy_control       - Run shared analysis with the healthy mask.
%   run_patient               - Run shared analysis with the patient mask.
%   run_hypoxia_bold_core     - Cohort-independent numerical analysis.
%
% Preprocessing
%   filter_bold_timeseries    - Frequency filtering and 12-volume rloess.
%   DelayDetermination_40_O2old - Initial voxelwise delay estimation.
%
% Normative analysis
%   compute_g_scalar          - Subject-level global response scalar.
%   fit_g_adjusted_model      - Corrected voxelwise OLS model.
%   apply_g_adjusted_model    - Expected, difference and z maps.
%   compute_voxelwise_icc     - Paired ICC(2,1) and ICC(3,1).
%
% Utilities
%   default_parameters        - Manuscript parameter values.
%   sensitivity_configurations - Seven-configuration sensitivity grid.
%   save_maps_mat             - Atomic MAT-file output.
%   write_maps_nifti_spm      - Primary NIfTI outputs using SPM12.
%
% Documentation
%   docs/INPUTS_AND_FILE_STRUCTURE.md - Required preprocessing and file tree.
