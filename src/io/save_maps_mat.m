% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function save_maps_mat(outputFile, maps, metadata)
%SAVE_MAPS_MAT Save computed maps and provenance metadata atomically.

if nargin < 3; metadata = struct(); end
outputDirectory = fileparts(outputFile);
if ~isempty(outputDirectory) && ~isfolder(outputDirectory); mkdir(outputDirectory); end
if isempty(outputDirectory); outputDirectory = pwd; end
temporaryFile = [tempname(outputDirectory), '.mat'];
save(temporaryFile, 'maps', 'metadata', '-v7.3');
movefile(temporaryFile, outputFile, 'f');
end
