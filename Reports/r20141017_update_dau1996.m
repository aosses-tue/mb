function r20141017_update_dau1996
% function r20141017_update_dau1996
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 15/10/2014
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 15/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

%% SDT Theory
% See Stanislaw1999, Tables 5 and 6

x = 1:6;
corr_tr = 0.5; % To avoid 1 or 0 H-FA rates
trial_sig = [0+corr_tr 4 8 8 12 18];
trial_noise = [8 6 1 3 7 0+corr_tr];
N_sig   = sum(trial_sig); % number of trials
N_noise = sum(trial_noise); % number of trials

H   = 1 - cumsum(trial_sig)/N_sig;
FA  = 1 - cumsum(trial_noise)/N_noise;

d = dprime(H,FA);

% Table 5:
disp('Stanislaw1999, Table 5: ')
var2latex([x' trial_sig' H' trial_noise' FA'])

% Table 6:
disp('Stanislaw1999, Table 6 (only d''): ')
var2latex([x' d'])

%% Adaptation loops:


%% SDT applied to 2 distributions
% To obtain the distributions:

close all

options.bPlot = 1;
options.dB_SPL = 85;
outs85 = demo_dau1996b(options);

options.bPlot = 0;
options.dB_SPL = 60;
outs60 = demo_dau1996b(options);

idx = outs85.out_stim1.idx;
t = outs85.out_stim1.t;

%%

% mue1 = optimaldetector(outs85.out_stim3.template_no_norm(:,idx),outs60.out_stim3.template_no_norm(:,idx))
% mue1 = optimaldetector(outs85.out_stim3.template(:,idx),outs60.out_stim3.template(:,idx))


%%
figure;
plot(   t,outs85.out_stim3.template(:,idx), ...
        t,outs60.out_stim3.template(:,idx) );

figure;
plot(   t,outs85.out_stim3.template_no_norm(:,idx), ...
        t,outs60.out_stim3.template_no_norm(:,idx) );
legend('85 dB SPL','60 dB SPL')
grid on

Ncentres = 1000;
[pdf1 centres1] = Probability_density_function(x1,Ncentres);
[pdf2 centres2] = Probability_density_function(x2,Ncentres);

x_pdf = 0:.01:400;
pdf1 = interp1(centres1,pdf1,x_pdf);
pdf2 = interp1(centres2,pdf2,x_pdf);

figure; 
plot(x_pdf,pdf1), hold on
plot(x_pdf,pdf2,'r');
grid on

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
