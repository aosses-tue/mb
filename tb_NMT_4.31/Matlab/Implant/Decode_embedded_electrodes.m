function [electrodes, modes] = Decode_embedded_electrodes(es, ms)

% Decode_embedded_electrodes: Decode embedded protocol electrode parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODE = Implant_modes;

electrodes = es;

is_stim_pulse = (es < 23);

mp1_present  = any(is_stim_pulse & (ms == 24));
mp2_present  = any(is_stim_pulse & (ms == 25));
mp12_present = any(is_stim_pulse & (ms == 30));

% Only one stimulation mode should occur:
if (mp1_present + mp2_present + mp12_present > 1)
	error('Mode must be constant')
end

if mp1_present

	modes = MODE.MP1;
	
elseif mp2_present

	modes = MODE.MP2;

else

	modes = MODE.MP1+2;

	mp12_idles = (es == 24) & (ms == 25);
	electrodes(mp12_idles) = 0;

end
