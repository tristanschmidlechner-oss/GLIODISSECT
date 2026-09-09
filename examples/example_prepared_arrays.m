% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
% Example using a prepared subject MAT file.
packageRoot = fileparts(fileparts(mfilename('fullpath')));
cd(packageRoot);
startup;

input = load('prepared_subject.mat');
if isfield(input, 'lesionMask')
    maps = run_patient(input);
else
    maps = run_healthy_control(input);
end

metadata = struct('matlab_version', version, 'parameters', default_parameters());
save_maps_mat('hypoxia_bold_maps.mat', maps, metadata);
