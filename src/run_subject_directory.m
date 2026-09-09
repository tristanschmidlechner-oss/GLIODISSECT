% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function maps = run_subject_directory(subjectDirectory, outputDirectory, cohort, filteringMode)
%RUN_SUBJECT_DIRECTORY Recompute primary EPI maps from explicitly selected inputs.
% Inputs may be read-only. outputDirectory must be a new directory outside
% subjectDirectory.
% filteringMode defaults to 'fresh'; 'historical_cache' explicitly uses the
% subject's root bBOLD.mat for filtering-independent core-equivalence tests.
% Cohort is 'healthy_control' or 'patient'. No ROI or plotting tail is run.
if nargin < 4, filteringMode = 'fresh'; end
assert(ismember(filteringMode,{'fresh','historical_cache'}), ...
    'HypoxiaBOLD:FilteringMode','Unknown filtering mode.');
subjectPath = char(java.io.File(char(subjectDirectory)).getCanonicalPath());
outputPath = char(java.io.File(char(outputDirectory)).getCanonicalPath());
assert(~strcmpi(subjectPath,outputPath) && ...
    ~startsWith([outputPath filesep],[subjectPath filesep],'IgnoreCase',ispc), ...
    'HypoxiaBOLD:UnsafeOutput','Output directory must be outside the input subject directory.');
assert(~isfolder(outputDirectory),'HypoxiaBOLD:ExistingOutput','Output directory already exists.');
assert(ismember(cohort,{'healthy_control','patient'}),'HypoxiaBOLD:Cohort','Unknown cohort.');
mkdir(outputDirectory);
boldFiles = dir(fullfile(subjectDirectory,'BOLD','32sr*.nii'));
assert(~isempty(boldFiles),'HypoxiaBOLD:MissingBOLD','No 32sr BOLD files.');
[~,order] = sort({boldFiles.name}); boldFiles=boldFiles(order);
paths = arrayfun(@(f) fullfile(f.folder,f.name),boldFiles,'UniformOutput',false);
headers = spm_vol(char(paths));
input = struct();
tissueNames = {'greyMatter','whiteMatter','CSF'};
for k=1:3
    tissueFile = exactly_one(fullfile(subjectDirectory,'T1',sprintf('reslicec%ds*.nii',k)));
    h = spm_vol(tissueFile);
    assert(isequal(h.dim,headers(1).dim) && max(abs(h.mat(:)-headers(1).mat(:)))<1e-5,...
        'HypoxiaBOLD:Geometry','Tissue and BOLD geometry differ; no registration guessed.');
    input.(tissueNames{k}) = spm_read_vols(h);
end
parameters = default_parameters();
if strcmp(cohort,'patient')
    input.lesionMask = prepare_patient_lesion(subjectDirectory,outputDirectory,headers(1));
    mask = build_patient_mask(input.greyMatter,input.whiteMatter,input.CSF,input.lesionMask,input,parameters);
else
    mask = build_healthy_mask(input.greyMatter,input.whiteMatter,input.CSF,input,parameters);
end
if strcmp(filteringMode,'historical_cache')
    cached = load(fullfile(subjectDirectory,'bBOLD.mat'),'bBOLD');
    input.bBOLD = cached.bBOLD; clear cached;
    assert(isequal(size(input.bBOLD),[headers(1).dim numel(headers)]), ...
        'HypoxiaBOLD:CacheGeometry','Cached BOLD dimensions differ from source images.');
    fprintf('Using historical cached bBOLD: %d scans, %d analysis voxels\n',size(input.bBOLD,4),nnz(mask));
else
    raw = spm_read_vols(headers);
    % The confirmed batches replace zero masked BOLD samples with NaN.
    raw(raw==0) = NaN;
    fprintf('Filtering %d scans, %d analysis voxels\n',size(raw,4),nnz(mask));
    input.bBOLD = filter_bold_timeseries(raw,mask,1.8,12);
    clear raw;
end
timing = prepare_oxygen_timing(exactly_one(fullfile(subjectDirectory,'Respiract files','*EndTidal*.txt')),...
    fullfile(subjectDirectory,'Respiract files','Events.xlsx'),input.bBOLD,1.8);
save(fullfile(outputDirectory,'timing.mat'),'timing','filteringMode','subjectDirectory');
fprintf('Computing voxel lag; initial global lag %g\n',timing.initialGlobalLag);
input.brainLag = DelayDetermination_40_O2old(timing.interpolatedO2,input.bBOLD);
input.durationRA = timing.durationRA;
input.initialGlobalLag = timing.initialGlobalLag;
input.rowScans = 1:size(input.bBOLD,4);
input.shiftedO2 = @(lag) resample_corrected_oxygen(timing,lag);
fprintf('Computing shared DTP/DTR/steady-state core\n');
if strcmp(cohort,'patient'), maps=run_patient(input); else, maps=run_healthy_control(input); end
write_maps_nifti_spm(maps,paths{1},outputDirectory);
save(fullfile(outputDirectory,'maps.mat'),'maps','-v7.3');
fprintf('Subject complete: %s\n',cohort);
end

function path = exactly_one(pattern)
files = dir(pattern);
assert(numel(files)==1,'HypoxiaBOLD:AmbiguousInput','Expected exactly one input matching %s; found %d.',pattern,numel(files));
path = fullfile(files.folder,files.name);
end

function lesion = prepare_patient_lesion(subjectDirectory,outputDirectory,boldHeader)
% savetumorvoxels provenance: NN-resliced labels 1/2/3 AND positive T1.
segmentation = fullfile(subjectDirectory,'segmentationoutput','oncohabitats','native','Segmentation.nii');
assert(isfile(segmentation),'HypoxiaBOLD:MissingSegmentation','Native segmentation is required.');
t1 = exactly_one(fullfile(subjectDirectory,'T1','reslices*.nii'));
h = spm_vol(t1);
assert(isequal(h.dim,boldHeader.dim) && max(abs(h.mat(:)-boldHeader.mat(:)))<1e-5,...
    'HypoxiaBOLD:Geometry','Resliced T1 and BOLD geometry differ.');
localSeg = fullfile(outputDirectory,'Segmentation.nii');
copyfile(segmentation,localSeg);
% spm_reslice is the write operation used by spm_run_coreg; write only the
% output-directory segmentation, with confirmed nearest-neighbour options.
flags = struct('interp',0,'wrap',[0 0 0],'mask',0,'which',1,'mean',0,'prefix','rB');
spm_reslice(char(t1,localSeg),flags);
labels = spm_read_vols(spm_vol(fullfile(outputDirectory,'rBSegmentation.nii')));
lesion = ismember(labels,[1 2 3]) & spm_read_vols(h)>0;
save(fullfile(outputDirectory,'segmentation_mask_evidence.mat'),'labels','lesion');
end
