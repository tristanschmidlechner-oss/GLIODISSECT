% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function test_historical_writer
%TEST_HISTORICAL_WRITER Verify the published map-specific encoding contract.
fixtureParent=fullfile(fileparts(fileparts(mfilename('fullpath'))),'test_outputs');
if ~isfolder(fixtureParent), mkdir(fixtureParent); end
root=tempname(fixtureParent); mkdir(root);
reference=struct('fname',fullfile(root,'reference.nii'),'dim',[2 2 2],...
    'dt',[16 0],'mat',eye(4),'pinfo',[1;0;0],'n',[1 1],'descrip','synthetic writer fixture');
spm_write_vol(reference,ones(2,2,2));
values=reshape([0 1 -1 NaN 2 3 4 5],[2 2 2]);
maps=struct('DTP_step1',values,'DTP_step2',values,'DTR_step1',values,...
    'DTR_step2',values,'corrected_delay',values,'steady_state_step1',values,...
    'steady_state_step2',values,'steady_state_O2norm_step1',values,'steady_state_O2norm_step2',values,...
    'stimulus_average_step1',values,'stimulus_average_step2',values,...
    'stimulus_average_O2norm_step1',values,'stimulus_average_O2norm_step2',values);
outputs=write_maps_nifti_spm(maps,reference.fname,fullfile(root,'outputs'));
assert(numel(outputs)==13,'All primary and comparator outputs must be written.');
for k=1:numel(outputs)
    got=spm_read_vols(spm_vol(outputs{k}));
    [~,name]=fileparts(outputs{k});
    if ismember(name,{'DTP_O2_map_step2','lag_i10_mask'}), expected=values;
    else, expected=values; expected(expected==0)=NaN; end
    assert(isequaln(got,expected),'Historical encoding mismatch for %s',name);
end
try
    write_maps_nifti_spm(maps,reference.fname,fullfile(root,'outputs'));
    error('HypoxiaBOLD:TestFailed','Existing outputs were not protected.');
catch err
    assert(strcmp(err.identifier,'HypoxiaBOLD:ExistingMap'));
end
fprintf('Historical writer and no-overwrite regression passed. Fixtures: %s\n',root);
end
