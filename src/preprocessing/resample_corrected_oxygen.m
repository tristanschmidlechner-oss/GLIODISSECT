% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function oxygen = resample_corrected_oxygen(timing, globalDelay)
%RESAMPLE_CORRECTED_OXYGEN Exact i10-aligned sampling from both batch scripts.
shift = globalDelay*timing.slopeAndIntercept(1);
ib = find(timing.time>timing.TimeOfVol(1)-shift,1,'first');
if isempty(ib), ib=1; end
ib = max(ib-1,1);
ie = find(timing.time>timing.TimeOfVol(end)-shift,1,'first');
if isempty(ie), ie=numel(timing.time); end
segment = timing.time(ib:ie);
grid = linspace(segment(1),segment(end),timing.numberOfScans)';
oxygen = interp1(segment,timing.oxygen(ib:ie),grid,'pchip','extrap')';
oxygen = oxygen(1:timing.numberOfScans);
end
