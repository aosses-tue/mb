function Hweight_lp = Get_Hweight_fluctuation2014(fs,outdir)
% function Hweight = Get_Hweight_fluctuation2014(fs,outdir)
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

if nargin == 0
    close all
    fs = 44100;
end

if nargin < 2
    outdir = 'D:\MATLAB_git\Psychoacoustics\FluctuationStrength_Garcia\private\';
end

% Design parameters of band-pass filter
sf1 = 0.5; % 0.5
pf1 = 3.1; % 2
pf2 = 12; % Hz % 8
sf2 = 20;
passAtt1 = 17.5;
passAtt2 = 14;

Hweight_lp = designfilt(   'lowpassiir', ...
                        'PassbandFrequency'  , pf2, ...
                        'StopbandFrequency'  , sf2, ...
                        'PassbandRipple'      , 3, ...
                        'StopbandAttenuation', passAtt2, ... % 100
                        'SampleRate'          , fs);
 
Hweight_hp = designfilt(   'highpassiir', ...
                        'StopbandFrequency'  , sf1, ...
                        'PassbandFrequency'  , pf1, ...
                        'StopbandAttenuation', passAtt1, ...
                        'PassbandRipple'      , 3, ...
                        'SampleRate'          , fs);

% HweightCACA = designfilt(   'bandpassiir', ...
%                         'StopbandFrequency1'  , sf1, ...
%                         'PassbandFrequency1'  , pf1, ...
%                         'PassbandFrequency2'  , pf2, ...
%                         'StopbandFrequency2'  , sf2, ...
%                         'StopbandAttenuation1', passAtt1, ...
%                         'PassbandRipple'      , 3, ...
%                         'StopbandAttenuation2', passAtt2, ... % 100
%                         'SampleRate'          , fs);

Hweight_HP = Hweight_hp.Coefficients;
Hweight_LP = Hweight_lp.Coefficients;
save([outdir 'Hweight-44100-Hz-LP.mat'],'Hweight_LP');
save([outdir 'Hweight-44100-Hz-HP.mat'],'Hweight_HP');

% figure;
% freqz(Hweight_lp,2^16);
% xlim([0.5 32]/(fs/2))
% 
% figure;
% freqz(Hweight_hp,2^16);
% xlim([0.5 32]/(fs/2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
