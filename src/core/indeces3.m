% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function [IndexBOLD90b,IndexBOLD10b, meanbaseline,....
    meanstep] = indeces3(Startsteptime, IndexBOLD10, IndexBOLD90, Endsteptime, fData,baselineshiftless, shift)
% Bas van Niftrik, Marco Piccirelli, Jorn Fierstra
% Department of Neurosurgery, University Hospital Zurich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The goal of this script is to determine the 10% and 90% index. It calculates
% the mean BOLD during the baseline and stepchange and index that does not 
% change anymore. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information:
% 1. Parameter explanation
% 2. Parameter reading
% 3. Create a new baseline row and a steptime row
% 4. Determine the new mean over the baseline and stepchange. 
% 5. Determine the new shifted 10% and 90%
% 6. Find the new 10% and 90% index value
% 7. Determine the new indices for both the 10% and 90%
% 8. Check to make sure the found 10% is the right 10% found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter explanation: 
% baselineshift                 Shift left over the baseline to correct for
%                               the last know 10% Index and used for the
%                               stepchange also for the left shift.
% baselineshiftless             Duration of baseline without any shift.
% baselinetiumeshift            Row of indices during the baseline with the
%                               shift subtracted.
% BOLD10_2                      BOLDbaseline + 10% value to determine the
%                               10% border
% BOLDh10_2                     Indices of BOLD from the 90% index to the 
%                               baselineshiftless.
% BOLD10per_2                   10% difference between baseline and step
% BOLD90_2                      BOLDbaseline + 90% value to determine the
%                               90% border
% BOLDh90_2                     Indices from 90% and all indices with BOLD
%                               values higher than the 90% border
% BOLD90per_2                   90% difference between baseline and step
% BOLDbaseline_2                BOLD values during the shifted baseline
% BOLDstep_2                    BOLD values during the shifted stepchange
% fData                         Squeezed functional MRi timeserie 2D
% IndexBOLD10b                  Precise Index of the 10% BOLD value
% IndexBOLD90b                  Precise Index of the 90% BOLD value
% Meanboldbaseline              Mean BOLD over the shifted baseline
% Meanboldstep                  Mean BOLD over the shifted stepchange
% steptime                      Row of indices during the
%                               stepchange(default size is 40)
% steptimerightshift            the right shift to the 90% index over the
%                               stepchange
% steptimeshift                 Row of indices during the stepchange with the
%                               two shifts subtracted.
%% Create a new baseline row and a steptime row
% Because you know that the shift at the end of the baseline to the 10% and
% form the start to the stepchange to the 90% you know the minimal left
% shift of the baseline and the minimal right shift of the stepchange.
% Furthermore, you know the left shift at the end of the stepchange

baselinetimeshifted = 1:1:IndexBOLD10;

% create a new array of the stepchange.
steper = IndexBOLD10 + 45; % Extend the candidate step window by five volumes.

steptimeshift = IndexBOLD90:1:steper;

%% Determine the new mean over the baseline and stepchange. 

Boldbaseline_2 = fData(baselinetimeshifted);
meanboldbaselineold = mean(Boldbaseline_2);

Boldstep_2 = fData(steptimeshift);
meanboldstepold = mean(Boldstep_2);


%% Determine the new shifted 10% and 90%

BOLD10per_2 = (meanboldstepold-meanboldbaselineold).*0.1;
BOLD10_2c = BOLD10per_2 + meanboldbaselineold;

BOLD90per_2 = (meanboldstepold-meanboldbaselineold).*0.9;
BOLD90_2c = BOLD90per_2 + meanboldbaselineold;

%% find the new 10% and 90% index value
% Divide between voxels with a positive and a negative signal. 10% can
% either be higher or lower than the baseline.
BOLDh90_2 = -1;
if BOLD90_2c >= meanboldbaselineold
    BOLDh90_2 = find(gt(fData(baselineshiftless:1:steper),BOLD90_2c));
    
elseif BOLD90_2c < meanboldbaselineold
    BOLDh90_2 = find(lt(fData(baselineshiftless:1:steper),BOLD90_2c));
end
%% Determine the new indices for both the 10% and 90%
if ~isempty(BOLDh90_2)
    IndexBOLD90b = BOLDh90_2(1)+ fix(baselineshiftless);
else
    BOLDh90_2(1) = IndexBOLD90;
    IndexBOLD90b = BOLDh90_2(1)+ fix(baselineshiftless);
end


BOLDh10_2 = -1;
if BOLD10_2c >= meanboldbaselineold
    fData2 = fData(IndexBOLD90b:-1:baselineshiftless);
    BOLDh10_2 =find(lt(fData2,BOLD10_2c));
elseif BOLD10_2c < meanboldbaselineold
    fData2 = fData(IndexBOLD90b:-1:baselineshiftless);
    BOLDh10_2 = find(gt(fData2,BOLD10_2c));
end
if isempty(BOLDh10_2)
    IndexBOLD10b = fix(baselineshiftless);
else
    IndexBOLD10b = IndexBOLD90b-BOLDh10_2(1);
end


%% create new means
bas = IndexBOLD10b+40;
if (bas-IndexBOLD90b)>= 0
    steptimeshift2 = IndexBOLD90b:1:bas;
else
    steptimeshift2 = (IndexBOLD90b-5):(IndexBOLD90b+5);
end

fData2 = fData(1:1:IndexBOLD10b);
meanbaseline = mean(fData2);
step_2 = fData(steptimeshift2);
meanstep = mean(step_2);
end
