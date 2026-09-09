% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
%RUN_SUBJECT_DIRECTORY_EXAMPLE Minimal directory-workflow example.
% Replace the four paths below. The output directory must not yet exist.

packageDir = '/path/to/GLIODISSECT';
spmDir = '/path/to/spm12';
subjectDir = '/path/to/preprocessed_subject';
outputDir = '/path/to/new_output_directory';

% Choose 'healthy_control' or 'patient'.
cohort = 'healthy_control';

restoredefaultpath
addpath(packageDir)
startup
addpath(spmDir)

run(fullfile(packageDir,'tests','run_all_tests.m'))
maps = run_subject_directory(subjectDir, outputDir, cohort);

fprintf('Finished. Results were written to:\n%s\n', outputDir);
