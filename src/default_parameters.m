% Copyright (C) 2015-2026 Tristan Schmidlechner, Vittorio Stumpo,
% Bas van Niftrik, and Jorn Fierstra
% ASTRAN Lab, Department of Neurosurgery, University Hospital Zurich
% SPDX-License-Identifier: BSD-3-Clause
%
function parameters = default_parameters()
%DEFAULT_PARAMETERS Analysis parameters used in the manuscript.

parameters.TR_seconds = 1.8;
parameters.temporal_window_volumes = 12;
parameters.temporal_window_seconds = 21.6;
parameters.low_pass_hz = 0.125;
parameters.high_pass_hz = 0;
parameters.dtp_lower_fraction = 0.10;
parameters.dtp_upper_fraction = 0.90;
parameters.spatial_smoothing_fwhm_mm = 6;
parameters.tissue_probability_threshold = 0.80;
parameters.normative_minimum_coverage = 0.80;
parameters.deviation_z_threshold = 2;
parameters.residual_sd_floor_percentile = 5;
end
