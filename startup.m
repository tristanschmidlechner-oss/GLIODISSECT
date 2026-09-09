% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function startup()
%STARTUP Add the publication package source directories to the MATLAB path.

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'src', 'core'));
addpath(fullfile(rootDir, 'src', 'preprocessing'));
addpath(fullfile(rootDir, 'src', 'masks'));
addpath(fullfile(rootDir, 'src', 'normative'));
addpath(fullfile(rootDir, 'src', 'io'));
end
