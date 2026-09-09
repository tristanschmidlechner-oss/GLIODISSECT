% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function report = compare_primary_precision(generatedDirectory, referenceDirectory)
%COMPARE_PRIMARY_PRECISION Whole-volume finite-pattern and precision gate.
% Timing must be exact. Continuous maps permit <= one float32 ULP at each
% compared value, accounting for rounding on a float32 write. All voxels are
% tested; no mask intersection or high-correlation exception hides differences.
mapping = {
'DTP_step1','DTP_O2_map_step1_epiMasked.nii','DTP_O2_map_step1_epiMasked.nii';
'DTP_step2','DTP_O2_map_step2_epiMasked.nii','DTP_O2_map_step2.nii';
'DTR_step1','DTR_O2_map_step1_epiMasked.nii','DTR_O2_map_step1_epiMasked.nii';
'DTR_step2','DTR_O2_map_step2_epiMasked.nii','DTR_O2_map_step2_epiMasked.nii';
'corrected_delay','lag_i10_mask.nii','lag_i10_mask.nii';
'HypoxiaSS_step1','HypoxiaStep1Map_SS_epiMasked.nii','HypoxiaStep1Map_SS_epiMasked.nii';
'HypoxiaSS_step2','HypoxiaStep2Map_SS_epiMasked.nii','HypoxiaStep2Map_SS_epiMasked.nii';
'HypoxiaSS_O2norm_step1','HypoxiaStep1Map_SS_O2norm_epiMasked.nii','HypoxiaStep1Map_SS_O2norm_epiMasked.nii';
'HypoxiaSS_O2norm_step2','HypoxiaStep2Map_SS_O2norm_epiMasked.nii','HypoxiaStep2Map_SS_O2norm_epiMasked.nii'};
rows=cell(9,12);
for k=1:9
    generatedPath=fullfile(generatedDirectory,mapping{k,2});
    if k==2 && isfile(fullfile(generatedDirectory,'DTP_O2_map_step2.nii'))
        generatedPath=fullfile(generatedDirectory,'DTP_O2_map_step2.nii');
    end
    ga=spm_vol(generatedPath);
    rb=spm_vol(fullfile(referenceDirectory,mapping{k,3}));
    assert(isequal(ga.dim,rb.dim) && max(abs(ga.mat(:)-rb.mat(:)))<1e-5, ...
        'HypoxiaBOLD:ComparisonGeometry','Comparison geometry differs.');
    a=spm_read_vols(ga); b=spm_read_vols(rb);
    valid=isfinite(a) & isfinite(b);
    pattern=nnz(xor(isfinite(a),isfinite(b)));
    nonfinite=nnz(~valid & ~((isnan(a)&isnan(b)) | (a==b)));
    delta=abs(a(valid)-b(valid));
    ulp=max(double(eps(single(abs(a(valid))))),double(eps(single(abs(b(valid))))));
    violations=nnz(delta>ulp);
    exactMismatch=nnz(delta~=0);
    if k<=5, pass=exactMismatch==0; rule='exact timing';
    else, pass=violations==0; rule='one float32 ULP'; end
    pass=pass && pattern==0 && nonfinite==0 && ~isempty(delta);
    rmse=NaN; maximum=NaN; maxULP=NaN;
    if ~isempty(delta), rmse=sqrt(mean(delta.^2)); maximum=max(delta); maxULP=max(delta./ulp); end
    rows(k,:)={mapping{k,1},nnz(isfinite(a)),nnz(isfinite(b)),pattern,exactMismatch,...
        violations,rmse,maximum,maxULP,nonfinite,rule,pass};
end
report=cell2table(rows,'VariableNames',{'map','generated_finite','reference_finite', ...
    'finite_pattern_mismatch','nonexact_finite','float32_ulp_violations','rmse', ...
    'maximum_absolute_difference','maximum_float32_ulp','nonfinite_pattern_mismatch','rule','pass'});
end
