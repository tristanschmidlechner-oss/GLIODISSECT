% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
packageRoot = fileparts(fileparts(mfilename('fullpath')));
cd(packageRoot);
startup;
addpath(fullfile(packageRoot,'tests'));

static_package_check;
test_masks;
test_normative_model;
test_icc;
test_filter_shape;
test_historical_writer;

disp('HypoxiaBOLD publication-package tests passed.');
