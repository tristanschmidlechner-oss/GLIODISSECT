% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function writtenFiles = write_maps_nifti_spm(maps, referenceNifti, outputDirectory)
%WRITE_MAPS_NIFTI_SPM Write primary maps using an SPM NIfTI header.

if exist('spm_vol','file') ~= 2 || exist('spm_write_vol','file') ~= 2
    error('HypoxiaBOLD:SPMUnavailable', 'SPM12 must be on the MATLAB path.');
end
if ~isfolder(outputDirectory); mkdir(outputDirectory); end
reference = spm_vol(referenceNifti);
reference = reference(1);
reference.dt = [16 0];

mapping = { ...
    'DTP_step1','DTP_O2_map_step1_epiMasked.nii'; ...
    'DTP_step2','DTP_O2_map_step2.nii'; ...
    'DTR_step1','DTR_O2_map_step1_epiMasked.nii'; ...
    'DTR_step2','DTR_O2_map_step2_epiMasked.nii'; ...
    'corrected_delay','lag_i10_mask.nii'; ...
    'steady_state_step1','HypoxiaStep1Map_SS_epiMasked.nii'; ...
    'steady_state_step2','HypoxiaStep2Map_SS_epiMasked.nii'; ...
    'steady_state_O2norm_step1','HypoxiaStep1Map_SS_O2norm_epiMasked.nii'; ...
    'steady_state_O2norm_step2','HypoxiaStep2Map_SS_O2norm_epiMasked.nii'; ...
    'stimulus_average_step1','HypoxiaAvg_step1_epiMasked.nii'; ...
    'stimulus_average_step2','HypoxiaAvg_step2_epiMasked.nii'; ...
    'stimulus_average_O2norm_step1','HypoxiaAvg_step1_O2norm_epiMasked.nii'; ...
    'stimulus_average_O2norm_step2','HypoxiaAvg_step2_O2norm_epiMasked.nii'};

writtenFiles = cell(size(mapping,1),1);
for index = 1:size(mapping,1)
    field = mapping{index,1};
    if ~isfield(maps, field); continue; end
    output = reference;
    output.fname = fullfile(outputDirectory, mapping{index,2});
    if isfield(output,'private') && isfield(output.private,'dat')
        output.private.dat.fname = output.fname;
    end
    data = double(maps.(field));
    % Both confirmed batches encode zero masked response/timing values as
    % NaN. DTP step 2 is written directly; corrected delay preserves zero
    % through the historical +0.001/-0.001 writer convention.
    if ~ismember(field, {'DTP_step2','corrected_delay'})
        data(data == 0) = NaN;
    end
    assert(~isfile(output.fname),'HypoxiaBOLD:ExistingMap', ...
        'Refusing to overwrite an existing map: %s',output.fname);
    spm_write_vol(output, data);
    writtenFiles{index} = output.fname;
end
writtenFiles = writtenFiles(~cellfun(@isempty,writtenFiles));
end
