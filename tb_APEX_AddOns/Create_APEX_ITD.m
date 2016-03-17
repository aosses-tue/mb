function y = Create_APEX_ITD
% function y = Create_APEX_ITD
%
% 1. Description:
%       Creates an adaptive procedure, using a similar structure to that in
%       gapdetection-mod (but updated to APEX v3.04
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 15/03/2016
% Last update on: 15/03/2016 
% Last use on   : 15/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

p = [];
% p.script='a3itd';
% 
% p=ef(p, 'experiment_file_prefix', 'a3itd');
p   = ef(p,'targetpath', [Get_TUe_paths('outputs') 'ITDac' delim]);
p   = ef(p,'stimdir', 'ITDac');
p.targetpath = makedirend(p.targetpath,1);
p.stimpath   = makedirend([p.targetpath p.stimdir],1);

Mkdir(p.targetpath);
if (~exist(p.stimpath))
    display('Maybe you are overwriting some information. Press any key to proceed...');
    pause();
end
Mkdir(p.stimpath);

p=ef(p,'fs', 44100);

p=ef(p,'itds', [-500:100:500]);  % us
p=ef(p,'ild' , 0);               % ILD to be introduced

% PROC_AFC=1;
% PROC_LR =2;
% PROC_TRAIN=3;
% 
% p=ef(p,'procedure', PROC_LR);          % procedure to be used
% p=ef(p,'add_reference', 0);             % present reference stimulus first
% 
p=ef(p,'sp_L', struct);                 % stimulus parameters Left (for genstimulus)
p=ef(p,'sp_R', struct);

p.sp_L=ef(p.sp_L,'stimulus_type', 1);
p.sp_R=ef(p.sp_R,'stimulus_type', 1);

p.sp_L=ef(p.sp_L,'freq', 2000);
p.sp_R=ef(p.sp_R,'freq', 2000);
 
% p.sp_L=ef(p.sp_L,'transpose', 0);
% p.sp_R=ef(p.sp_R,'transpose', 0);
 
% p.sp_L=ef(p.sp_L,'noise_bandwidth', 100);
% p.sp_R=ef(p.sp_R,'noise_bandwidth', 100);
 
p.sp_L=ef(p.sp_L,'stimulus_len',0.5);
p.sp_R=ef(p.sp_R,'stimulus_len',0.5);
 
% datablocks  = [];
% stimuli     = [];
% trials      = [];
% connections = [];

lvl = 65;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Creates calibration tone:
% cal = p.sp_L;
% cal.stimulus_len=10; % make long masker for calibration
% cal.ramp_type=0;
% cal.ramp=0;
% calibfilename=[ num2str((p.sp_L.modfreq)) '_calib.wav'];
% % calib = genstimulus(cal);
% calib = Create_sin(cal.freq,cal.stimulus_len,p.fs);

% calib = setdbspl(calib,lvl);
% Wavwrite(calib,p.fs, [p.targetpath calibfilename]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if (p.procedure==PROC_LR)
%     % create silence datablocks
%     datablocks=[datablocks a3datablock('silence', 'silence:300', 'soundcard',2)];
% end

for itd=p.itds
    l = genstimulus_getname(p.sp_L); 
    r = genstimulus_getname(p.sp_R);
    filename=[l '_' r '_itd' num2str(itd) '.wav'];
    totalfilename=[p.stimpath filename];

    if (itd<0)
        p.sp_R.shift=0;
        p.sp_L.shift=-itd/1e6;
    else
        p.sp_R.shift=itd/1e6;
        p.sp_L.shift=0;
    end
    
    left  = Create_sin(p.sp_L.freq, p.sp_L.stimulus_len, p.fs); % genstimulus(p.sp_L);
    right = Create_sin(p.sp_R.freq, p.sp_R.stimulus_len, p.fs); % genstimulus(p.sp_R);
    left  = setdbspl(left,lvl);
    right = setdbspl(right,lvl);
        
    left  = il_shift(left ,p.sp_L.shift,p.fs);
    right = il_shift(right,p.sp_R.shift,p.fs);     
        
    d=length(right)-length(left);
    if (d>0)
        left=[left; zeros(d,1)];
    elseif (d<0)
        right=[right; zeros(-d,1)];
    end

    % introduce ILD if necessary
    if (p.ild)
        right=right*10^(p.ild/20);
    end

    data=[left right];

    Wavwrite(data, p.fs, totalfilename);
    % fixparam.itd=itd;
   
end

disp('')

if bDiary
	diary off
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=il_itdstring(str, itd)
        if (itd==-0), itd=0; end;
    result=sprintf('%s%g', str,itd);
    return;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outsig = il_shift(insig,shift,fs)

nsamp=round(shift*fs);
if (abs(shift/fs)>1e6)     % microsecond accuracy
    warning(['p.shift rounded to nearest integer number of samples: ' num2str(nsamp) ]);
end
outsig=[zeros(nsamp,1); insig];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
