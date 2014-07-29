function AMT_running_demos
% function AMT_running_demos
%
%   Demonstration of model stages
%     DEMO_ABSOLUTETHRESHOLD - Absolute thresholds of hearing
%     DEMO_ADAPTLOOP         - Adaptation loops.
%     DEMO_DRNL              - Widening of filters in the DRNL
%
%   Demonstration of full models
%     DEMO_LINDEMANN1986     - Lindemann binaural model.
%     DEMO_BAUMGARTNER2013   - Sagittal-plane sound localization (Baumgartner et al., 2013)
%     DEMO_HOHMANN2007       - Hohmann 2007 filterbank
%     DEMO_TAKANEN2013       - Takanen et al. (2013) model
%     DEMO_ZILANY2013        - Zilaney model
%     DEMO_VERHULST2012      - Verhulst et al. (2012)
%
%  1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
% % Before running (Unix-based):
%   misc = Get_TUe_paths([], 1);
%   ltfatstart
%   ltfatmex
%
% % Before running (Win-based):
%   Start_TUe;
%   AMT_running_demos
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: March 2014
% Last update: 14/05/2014 % Update this date manually
% Last used: 14/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.bSave = 0;
info.outputs = Get_TUe_paths('outputs');

Mkdir(info.outputs); % Creates directory only if it does not exist

close all

count_figs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Absolute threshold

h = demo_absolutethreshold;

% Figure 1 (1): Thresholds of hearing by standard
% Figure 2 (2): Absolute thesholds for the ER2A and the Sennheiser HDA-200 
%           earphones are provided up to 16 kHz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. demo_adaptloop: non-linear adaptation

% Dudas:
%       What is 'overshoot limiting'

h(end+1:end+2) = demo_adaptloop;

%   Figure 1 (3): Clean test signal
%   Figure 2 (4): Noisy test signal

% Ref 25: J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing model based on contralateral inhibition. I. Model structure. J. Acoust. Soc. Am., 110:1074-1088, August 2001.
% Ref 35: T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the effective signal processing in the auditory system. I. Model structure. J. Acoust. Soc. Am., 99(6):3615-3622, 1996a.
% D. Pueschel. Prinzipien der zeitlichen Analyse beim Hoeren. PhD thesis, Universitaet Goettingen, 1988. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. demo_drnl: DRNL models the BM non-linearity. DRNL is followed by inner
%               cells envelope extraction (half-wave rectification + low-pass filtering)

[h(end+1:end+3), output] = demo_drnl;

if info.bSave % to generate wav-files
     Wavwrite(output.signal_at_lvl1,output.Fs,[info.outputs 'signal_at_lvl1'])
     Wavwrite(output.signal_at_lvl2,output.Fs,[info.outputs 'signal_at_lvl2'])
     Wavwrite(output.signal_at_lvl3,output.Fs,[info.outputs 'signal_at_lvl3'])
end

% Figure 1 (5), 2 (6), 3 (7): Input at different levels:
%   - Plots with filterbank of about 300 band-resolution
%   - Not clear the amplitude units

% E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank. JASA. 110:3107-3118, 2001.

%% Full models:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. LINDEMANN1986
%   type: binaural model
h(end+1:end+2) = demo_lindemann1986;

% 5. Lindemann's model with 'own data': see Run_HRIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. DEMO_VERHULST2012

% demo_verhulst2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if info.bSave
    for i = 1:length(h)

        Saveas(h(i),[info.outputs 'AMT-fig-' num2str(count_figs)]); 
        count_figs = count_figs+1; 

    end
end

end
% EOF