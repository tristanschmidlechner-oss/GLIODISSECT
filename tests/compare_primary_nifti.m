% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function metrics = compare_primary_nifti(generatedDirectory,referenceDirectory,mask)
%COMPARE_PRIMARY_NIFTI Compare all nine primary maps without resampling.
% Report finite-domain differences AND missing-pattern differences; finite
% intersection correlations alone cannot establish numerical reproduction.
pairs = {
'DTP_step1','DTP_O2_map_step1_epiMasked.nii','DTP_O2_map_step1_epiMasked.nii';
'DTP_step2','DTP_O2_map_step2_epiMasked.nii','DTP_O2_map_step2.nii';
'DTR_step1','DTR_O2_map_step1_epiMasked.nii','DTR_O2_map_step1_epiMasked.nii';
'DTR_step2','DTR_O2_map_step2_epiMasked.nii','DTR_O2_map_step2_epiMasked.nii';
'corrected_delay','lag_i10_mask.nii','lag_i10_mask.nii';
'HypoxiaSS_step1','HypoxiaStep1Map_SS_epiMasked.nii','HypoxiaStep1Map_SS_epiMasked.nii';
'HypoxiaSS_step2','HypoxiaStep2Map_SS_epiMasked.nii','HypoxiaStep2Map_SS_epiMasked.nii';
'HypoxiaSS_O2norm_step1','HypoxiaStep1Map_SS_O2norm_epiMasked.nii','HypoxiaStep1Map_SS_O2norm_epiMasked.nii';
'HypoxiaSS_O2norm_step2','HypoxiaStep2Map_SS_O2norm_epiMasked.nii','HypoxiaStep2Map_SS_O2norm_epiMasked.nii'};
rows=cell(size(pairs,1),12);
for k=1:size(pairs,1)
    generatedPath=fullfile(generatedDirectory,pairs{k,2});
    if k==2 && isfile(fullfile(generatedDirectory,'DTP_O2_map_step2.nii'))
        generatedPath=fullfile(generatedDirectory,'DTP_O2_map_step2.nii');
    end
    ga=spm_vol(generatedPath);
    rb=spm_vol(fullfile(referenceDirectory,pairs{k,3}));
    assert(isequal(ga.dim,rb.dim) && max(abs(ga.mat(:)-rb.mat(:)))<1e-5,...
        'HypoxiaBOLD:ComparisonGeometry','Map geometry differs; no comparison resampling allowed.');
    a=spm_read_vols(ga); b=spm_read_vols(rb);
    common=mask & isfinite(a) & isfinite(b);
    missing=nnz(mask & xor(isfinite(a),isfinite(b)));
    d=a(common)-b(common);
    r=NaN; rmse=NaN; med=NaN; maximum=NaN; identical=NaN; within=NaN;
    if ~isempty(d)
        rmse=sqrt(mean(d.^2)); med=median(abs(d)); maximum=max(abs(d));
        identical=100*mean(d==0); within=100*mean(abs(d)<=1);
        if numel(d)>1, c=corrcoef(a(common),b(common)); r=c(1,2); end
    end
    rows(k,:)={pairs{k,1},nnz(mask & isfinite(a)),nnz(mask & isfinite(b)),nnz(common),missing,...
        r,rmse,med,maximum,identical,within,~isempty(d) && all(d==0) && missing==0};
end
metrics=cell2table(rows,'VariableNames',{'map','generated_finite','reference_finite','common_finite',...
    'finite_pattern_mismatch','pearson_r','rmse','median_absolute_difference','maximum_absolute_difference',...
    'percent_exact','percent_within_one','exact_match'});
end
