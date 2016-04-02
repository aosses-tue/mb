function makeitd_ac(stage)
% function makeitd_ac(stage)
%
% 1. Description:
%       Run this script having it in the current folder
%       p.procedure = 2 % Constant stimuli
% 
% 2. Stand-alone example:
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    stage = 1;
end

p.sp_L.normalize_rms= -20;
p.sp_R.normalize_rms= -20;
p.sp_L.len          =   0.1; % s
p.sp_R.len          =   0.1; % s
p.procedure         =   2;
p.ild               =  -12; % -4 dB
p.sp_L.freq         = 2000; % Hz
p.sp_R.freq         = 2000; % Hz

p.targetpath = 'stimuli'; %'/Anneke/ITDac/stimuli/';

if p.ild ~= 0
    p.targetpath = sprintf('%s-ILD-%.0f-dB',p.targetpath,p.ild);
end

p.dirstimuli = p.targetpath;

Mkdir(p.targetpath)

type = 1;
if (type==1)            % transposed
    if(stage==1)  
        p.itds=[0:150:900];
        
    elseif(stage==2)
        p.calibration_profile='ITD_ac';
        p.itds=[-500 :100: 500];
        p.sp_L.stimulus_type=1;
        p.sp_L.freq=2400;
        p.sp_L.modfreq=400;
        p.sp_L.transpose=1;
        p.sp_L.ramp = 1;
        p.sp_R=p.sp_L;

        p.targetpath=[targetpath 'trans_ac' num2str(p.sp_L.freq) '_mod' num2str(p.sp_L.modfreq) '_stage' num2str(stage)];

    elseif(stage==10)
        p.calibration_profile='ITD_ac';
        p.itds=[-250 :50: 250];
        p.sp_L.stimulus_type=1;
        p.sp_L.freq=2400;
        p.sp_L.modfreq=400;
        p.sp_L.transpose=1;
        p.sp_L.ramp = 1;
        p.sp_R=p.sp_L;

        p.targetpath=[targetpath 'trans_ac' num2str(p.sp_L.freq) '_mod' num2str(p.sp_L.modfreq) '_stage' num2str(stage)];

    elseif(stage==11)
        p.calibration_profile='ITD_ac';
        p.itds=[-250 :50: 250];
        p.sp_L.stimulus_type=1;
        p.sp_L.freq=2400;
        p.sp_L.modfreq=600;
        p.sp_L.transpose=1;
        p.sp_L.ramp = 1;
        p.sp_R=p.sp_L;

        p.targetpath=[targetpath 'trans_ac' num2str(p.sp_L.freq) '_mod' num2str(p.sp_L.modfreq) '_stage' num2str(stage)];

    elseif(stage==20)
        p.calibration_profile='ITD_ac';
        p.itds=[-250 :50: 250];
        p.sp_L.stimulus_type=1;
        p.sp_L.freq=1200;
        p.sp_L.modfreq=400;
        p.sp_L.transpose=1;
        p.sp_L.ramp = 1;
        p.sp_R=p.sp_L;

        p.targetpath=[targetpath 'trans_ac' num2str(p.sp_L.freq) '_mod' num2str(p.sp_L.modfreq) '_stage' num2str(stage)];

    elseif(stage==21)
        p.calibration_profile='ITD_ac';
        p.itds=[-250 :50: 250];
        p.sp_L.stimulus_type=1;
        p.sp_L.freq=1200;
        p.sp_L.modfreq=600;
        p.sp_L.transpose=1;
        p.sp_L.ramp = 1;
        p.sp_R=p.sp_L;

        p.targetpath=[targetpath 'trans_ac' num2str(p.sp_L.freq) '_mod' num2str(p.sp_L.modfreq) '_stage' num2str(stage)];

    else
    error('Invalid type');
    end
end
a3itd_ac(p,0);

% p.procedure=3;  % training
% a3itd_ac(p);