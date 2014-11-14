function gain_dB = SPrint_sensitivity_gain_dB(sens)
% function gain_dB = SPrint_sensitivity_gain_dB(sens)
%
% Convert user sensitivity value on LCD to gain value in dB.
% Each step is 1.5 dB.
% Sensitivity value 0 is treated specially, and gives 0 dB gain.
%
% gain_dB = SPrint_sensitivity_gain_dB(sens)
%
% sens:     Use sensitivity value on LCD (range 0 - 20), default = 12
% gain_dB:  Gain in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SG_BASE = 22 * 0.375;
SG_STEP =  4 * 0.375;

gain_dB = SG_BASE + SG_STEP * sens;
gain_dB(sens == 0) = 0;
