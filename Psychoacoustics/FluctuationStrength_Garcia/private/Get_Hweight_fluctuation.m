function Hweight = Get_Hweight_fluctuation(fs)
% function Hweight = Get_Hweight_fluctuation(fs)
% 
% Returns the Hweight filter.
% 
% Inputs:
% params: Struct specifying filter characteristics.
% fs: Sampling frequency.
% 
% Outputs:
% Hweight: The digital filter.
% 
% Author: Rodrigo Garcia/Alejandro Osses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Design parameters of band-pass filter
sf1 = 0.5;
pf1 = 2;
pf2 = 8;
sf2 = 32;

try
	Hweight = designfilt(...
        'bandpassiir', ...
        'StopbandFrequency1', sf1, ...
        'PassbandFrequency1', pf1, ...
        'PassbandFrequency2', pf2, ...
        'StopbandFrequency2', sf2, ...
        'StopbandAttenuation1', 100, ...
        'PassbandRipple', 3, ...
        'StopbandAttenuation2', 100, ...
        'SampleRate', fs);
catch
    disp('Make sure you have a MATLAB version R2014b or above...')
    if fs ~= 44100
        error('Run this at 44100 Hz or otherwise use MATLAB R2014b or above');
    end
    load('Hweight-44100-Hz.mat');
    % Hweight = 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
