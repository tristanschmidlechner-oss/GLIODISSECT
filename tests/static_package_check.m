% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function static_package_check
%STATIC_PACKAGE_CHECK Verify function resolution after startup.

packageRoot = fileparts(fileparts(mfilename('fullpath')));
required = {'run_healthy_control','run_patient','run_hypoxia_bold_core', ...
    'run_subject_directory','prepare_oxygen_timing','resample_corrected_oxygen','default_parameters', ...
    'computeDTP','computeDTR','complrange10_90c','indeces3', ...
    'OxygenDoubleStep_SS','compute_stimulus_average', ...
    'filter_bold_timeseries','low_highPassFilter_g', ...
    'DelayDetermination_40_O2old','build_healthy_mask','build_patient_mask', ...
    'compute_g_scalar','fit_g_adjusted_model','apply_g_adjusted_model', ...
    'compute_voxelwise_icc','sensitivity_configurations','save_maps_mat', ...
    'write_maps_nifti_spm'};

for index = 1:numel(required)
    resolved = which(required{index}, '-all');
    if isempty(resolved)
        error('HypoxiaBOLD:MissingFunction', 'Cannot resolve %s.', required{index});
    end
    if ischar(resolved); resolved = cellstr(resolved); end
    packageMatches = startsWith(resolved, [packageRoot filesep]);
    if numel(resolved) ~= 1 || ~all(packageMatches)
        error('HypoxiaBOLD:AmbiguousFunction', ...
            '%s must resolve exactly once from this package; found %d locations.', required{index}, numel(resolved));
    end
end

core=fileread(fullfile(packageRoot,'src','run_hypoxia_bold_core.m'));
assert(nargin('run_hypoxia_bold_core')==7);
assert(~contains(lower(core),'fallback'));
assert(~contains(core,'median(validDelay)'));
assert(isempty(which('highpassfiltering_VS','-all')));
assert(isempty(which('low_highPassFilter_g_VS','-all')));
disp('Static function-resolution and shared-core contract checks passed.');
end
