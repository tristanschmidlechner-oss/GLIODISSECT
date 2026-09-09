% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function brainlagDiffO2 = DelayDetermination_40_O2old(InterpolatedO2, bBOLD)
%DELAYDETERMINATION_40_O2OLD Estimate voxelwise BOLD-to-O2 lag.
% By Bas van Niftrik and Marco Piccirelli; adapted for O2 by Vittorio Stumpo.
%
% Input:
%   Interpolated O2: O2 interpolated data _ with the extra time ... 
%           a vector with the same number of elements scans PLUS a time
%           bigger than the maximum expected delays. 
%   BOLD: A string the the BOLD data 4D-hdr-file name including path
%
% Output: 
%   brainlagdiffO2: A 3D matrix of delay for each voxel of the mask

%% Interpolated O2 data
O2 = (InterpolatedO2 - mean(InterpolatedO2)) ./ std(InterpolatedO2);

%% Read the BOLD data
MtxSizeX = size(bBOLD,1);
MtxSizeY = size(bBOLD,2);
MtxSizeZ = size(bBOLD,3);

%% Apply the mask if needed

%% For each voxel find the optimal delay for maximum correlation

brainlagDiffO2 = zeros(MtxSizeX,MtxSizeY,MtxSizeZ);
for ix = 1:MtxSizeX
    for iy = 1:MtxSizeY
        for iz = 1:MtxSizeZ
            fData2 = squeeze(bBOLD(ix,iy,iz,:));
            % Retain historical xcorr/max behavior for every voxel, including
            % invalid/constant series. The original zero-lag background feeds
            % the unmasked DTP-step-2 map.
            fData = detrend(fData2); 
            fData = (fData-mean(fData))./std(fData);
            [max1, imax1] = xcorr(fData,O2,'none');
            [~,I] = max(abs(max1((imax1(end)):(imax1(end)+15)))); %after end + 30
            brainlagDiffO2(ix,iy,iz) = I-1;
        end
    end
end

end
