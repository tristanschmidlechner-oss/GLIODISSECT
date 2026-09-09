% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [ output ] = low_highPassFilter_g( TimeEvol, input, FilterShape)
% out = low_highfilter(TimeEvol, [150 10 0.001 1000], 'butter', 21 )
% input = [150 10 0.001 1000];
% FilterShape = 'butter';
% FigIndex = 21;
%% input arguments
% TimeEvol is the voxel Bold Signal
lowPassFreq = input(1);
highPassFreq = input(2);
TR = input(3); % =T % Sample time, TR [s]
L = input(4); % Length of signal in samples (volumes, dynamics)
FilterChoice = FilterShape;


Fs = 1/TR; % Sampling frequency [Hz]
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
NFFT = L;




%% Filterdesign
switch FilterChoice
    case 'box'
        % Box Filter
        Filter = zeros(1,NFFT);
        NLowStopFreqSample = floor(highPassFreq/(Fs/L))+1;
        NHighStopFreqSample = ceil(lowPassFreq/(Fs/L))-1;
        Filter(NLowStopFreqSample:NHighStopFreqSample)=1;

    case 'smooth box'
        % Box Filter
        Filter = zeros(1,NFFT);
        NLowStopFreqSample = floor(highPassFreq/(Fs/L))+1;
        %NLowStopFreqSample = floor(highPassFreq/(Fs/L));
        NHighStopFreqSample = ceil(lowPassFreq/(Fs/L))-1;
        Filter(NLowStopFreqSample:NHighStopFreqSample)=1;
        RampFreq=lowPassFreq/2;
        %        FreqRange=min(floor(L*RampFreq/Fs), floor(L*2*highPassFreq/Fs))-1;
        FreqRange=max(5,floor(L*RampFreq/Fs));
        %         Smoother=gausswin(FreqRange,2.944)/norm(gausswin(FreqRange,2.944)); % relative weight: 1, 0.500, 0.062
        Smoother=gausswin(FreqRange,2.944*FreqRange*(5-1)/(5*(FreqRange-1)))/norm(gausswin(FreqRange,2.944*FreqRange*(5-1)/(5*(FreqRange-1))));
        FilterTemp = conv(Filter,Smoother,'same');
        Filter(floor(mean([NLowStopFreqSample,NHighStopFreqSample])):end)=FilterTemp(floor(mean([NLowStopFreqSample,NHighStopFreqSample])):end)/...
            FilterTemp(floor(mean([NLowStopFreqSample,NHighStopFreqSample]))+1);
    otherwise
        error('HypoxiaBOLD:UnsupportedFilter', 'Filter shape "%s" is not implemented.', FilterChoice)
end

%% Filtering

Y = fft( TimeEvol ,NFFT)/L;
FilteredY = Y.*Filter';

FilteredBoldSignal=ifft(FilteredY, 'symmetric');




%% outputing
output = L.*FilteredBoldSignal;
