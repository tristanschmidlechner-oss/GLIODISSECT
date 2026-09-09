% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function timing = prepare_oxygen_timing(endTidalFile, eventsFile, bBOLD, TR)
%PREPARE_OXYGEN_TIMING Retain the confirmed Gen4 protocol-9 timing equations.
% Derived from parametersreadingZurich_O2_new_batch, marfLinearfit_gen4O2_old,
% sampleBOLDO2_gen4_new, and durationO2_deoxy3. Plotting is omitted.
% See docs/PROVENANCE.md for source lineage and release scope.
N = size(bBOLD,4);
XX = textread(endTidalFile,'%f'); %#ok<DTXTRD>
assert(mod(numel(XX),23)==0,'HypoxiaBOLD:EndTidalFormat','Expected 23-field records.');
oxygen = XX(6:23:end);
time = XX(1:23:end)/1000;
Events = xlsread(eventsFile); %#ok<XLSRD>
Events(:,1) = Events(:,1)/1000;
[time,indices] = unique(time);
oxygen = oxygen(indices);
eventRows = 2:(size(Events,1)-1);
fit = polyfit(Events(eventRows,2)',Events(eventRows,1)',1);
TimeOfVol = fit(1)*[1 N]+fit(2);
ib = find(time>TimeOfVol(1),1)-1;
ie = find(time>TimeOfVol(end),1);
assert(~isempty(ib) && ib>=1 && ~isempty(ie),'HypoxiaBOLD:TimeCoverage','Recording does not cover the BOLD interval.');
grid = linspace(time(ib),time(ie),N)';
interpolated = interp1(time(ib:ie),oxygen(ib:ie),grid,'pchip','extrap')';
v = nanmean(nanmean(nanmean(bBOLD)));
v = reshape(v,1,[]);
brain = detrend(v)+mean(v);
N2 = (brain-mean(brain))/std(brain);
O2E = (interpolated-mean(interpolated))/std(interpolated);
[c,lags] = xcorr(N2,O2E,'none');
[~,II] = max(abs(c(lags(end):(lags(end)+15))));
lag = II-1;
first = 60/TR-fix((TimeOfVol(1)-Events(1))/fit(1))+lag+3;
duration = [first 60/TR 40/TR 60/TR N-(first+60/TR+40/TR+60/TR)];
timing = struct('initialGlobalLag',lag,'durationRA',duration,...
    'interpolatedO2',interpolated,'time',time,'oxygen',oxygen,...
    'TimeOfVol',TimeOfVol,'slopeAndIntercept',fit,'numberOfScans',N);
end
