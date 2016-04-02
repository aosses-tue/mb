function a3itd_ac(p,doskip)
% function a3itd_ac(p,doskip)
%
% 1. Description:
%       if doskip is 1, already existing stimuli will be kept
%
% 2. Stand-alone example:
%       makeitd_ac.m
% 
% Some dependencies:
%   genstimulus.m
%   a3connection.m
%   a3soundlevelmeter.m
%   a3corrector.m
%   a3amplifier.m
%   bool2xml.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
    doskip=1;
end

p.script='a3itd';

p=ef(p, 'experiment_file_prefix', 'a3itd');
p=ef(p,'targetpath', [Get_TUe_paths('outputs') 'ITDac\']);
Mkdir(p.targetpath)

p=ef(p,'dirstimuli', 'stimuli');

p=ef(p,'fs', 44100);

p=ef(p,'itds', [-500:100:500]);         % us
p=ef(p,'ild',   0);                     % ILD to be introduced

PROC_AFC=1;
PROC_LR=2;
PROC_TRAIN=3;

p=ef(p,'procedure', PROC_LR); % procedure to be used
p=ef(p,'add_reference', 0);   % present reference stimulus first

p=ef(p,'sp_L', struct);       % stimulus parameters Left (for genstimulus)
p=ef(p,'sp_R', struct);

p.sp_L=ef(p.sp_L,'stimulus_type', 1);
p.sp_R=ef(p.sp_R,'stimulus_type', 1);

p.sp_L=ef(p.sp_L,'len',100e-3);
p.sp_R=ef(p.sp_R,'len',100e-3);

p.sp_L=ef(p.sp_L,'freq', 3900);
p.sp_R=ef(p.sp_R,'freq', 3900);

p.sp_L=ef(p.sp_L,'modfreq', 300);
p.sp_R=ef(p.sp_R,'modfreq', 300);

p.sp_L=ef(p.sp_L,'transpose', 0);
p.sp_R=ef(p.sp_R,'transpose', 0);

p.sp_L=ef(p.sp_L,'noise_bandwidth', 100);
p.sp_R=ef(p.sp_R,'noise_bandwidth', 100);

p.targetpath=makedirend(p.targetpath,1);
p.stimpath=makedirend([p.targetpath],1);

datablocks=[];
stimuli=[];
trials=[];
connections=[];

% Calibration:
cal=p.sp_L;
cal.stimulus_len=10;        % make long masker for calibration
cal.ramp_type=0;
cal.ramp=0;
calibfilename=[ num2str((p.sp_L.modfreq)) '_calib.wav'];
% calib = genstimulus(cal);
calib = Create_sin(cal.freq,cal.stimulus_len,p.fs);

Wavwrite(calib,p.fs, [p.targetpath calibfilename]);

if (p.procedure==PROC_LR || p.procedure==PROC_TRAIN)
    % ensure all + and - itd's are available
    for c=1:length(p.itds)
        if (~length(find(p.itds==-p.itds(c))))
            p.itds=[p.itds -p.itds(c)];
        end
    end
end

datablocks = ['<uri_prefix>' p.dirstimuli '</uri_prefix>' lf];
if (p.procedure==PROC_LR)
    % create silence datablocks
    datablocks=[datablocks a3datablock('silence', 'silence:300', 'soundcard',2)];
end

if (~length(p.itds))
    error('No itd''s given');
end

for itd=p.itds
    l=genstimulus_getname(p.sp_L);
    % data = Create_sin(p.sp_L.freq,cal.stimulus_len,p.fs);
    
    r=genstimulus_getname(p.sp_R);
    filename=[l '_' r '_itd' num2str(itd) '_modF' num2str(p.sp_L.modfreq)  '.wav'];
    totalfilename=[p.stimpath filename];

    if (~exist(totalfilename) || ~doskip)
        if (itd>0) % ITD > 0, Right leads
            p.sp_R.shift=0;
            p.sp_L.shift=abs(itd)/1e6;
        else
            p.sp_R.shift=abs(itd)/1e6;
            p.sp_L.shift=0;
        end
            
    else
        error('Method not implemented');
    end

    if (p.sp_L.stimulus_type==4 && p.sp_R.stimulus_type==4 & 0)        %%FIXME
        nbs=p.sp_L;
        %nbs.

        nbs.shift_fs=1;
        nbs.shift_env=1;
        nbs.itd=itd;
        %nbs.noise_cutoff=p.sp_L.freq*2^(1/p.sp_L.noise_bandwidth/2)-p.sp_L.freq/2^(1/p.sp_L.noise_bandwidth/2);
        %nbs.noise_cutoff=nbs.noise_cutoff/2;
        nbs.noise_cutoff=nbs.noise_bandwidth/2;
        nbs.freq=p.sp_L.freq;

        [left,right]=gennoiseband_special(nbs);
    else

        left = Create_sin(p.sp_L.freq,p.sp_L.len,p.fs);% genstimulus(p.sp_L);
        right= Create_sin(p.sp_R.freq,p.sp_R.len,p.fs);% genstimulus(p.sp_R);

        left  = il_shift(left ,p.sp_L.shift,p.fs);
        right = il_shift(right,p.sp_R.shift,p.fs);
        
        d=length(right)-length(left);
        if (d>0)
            left=[left; zeros(d,1)];
        elseif (d<0)
            right=[right; zeros(-d,1)];
        end
    end

    % introduce ILD if necessary
    if (p.ild)
        right=right*10^(p.ild/20);
    end

    data=[left right];

    Wavwrite(data, p.fs, totalfilename);
% fixparam.itd=itd;

datablocks=[datablocks a3datablock(itdstring('datablock', itd), filename, 'soundcard') ];

if (p.procedure==PROC_LR )
    if (itd==0)
%          continue;
    end

    if (p.add_reference)
        dbs={itdstring('datablock', itd) 'silence' itdstring('datablock', -itd)};
    else
        dbs={itdstring('datablock', -itd)};
    end
    if (itd<0)
        answ='links';
    else
        answ='rechts';
    end
    trials=[trials a3trial( itdstring('trial', itd), 'screen', itdstring('stimulus', itd), answ)];

elseif (p.procedure==PROC_TRAIN)
    dbs={itdstring('datablock', -itd)};
    answ=sprintf('button%3.2f', itd);
    trials=[trials a3trial( itdstring('trial', itd), 'screen', itdstring('stimulus', itd), answ)];
else
    error('Unknown procedure type');
end
fixparam.itd=itd;
stimuli=[stimuli a3stimulus(itdstring('stimulus', itd), dbs, fixparam,struct,0) ];

connections=[connections a3connection(itdstring('datablock',itd),1,'soundcard', 1) lf ];
end

%% Calibration datablock
datablocks=[datablocks a3datablock(['datablock_calib'],calibfilename ,'soundcard')];
%%%%%%%%%%%%%%%%%%%%%%% calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%
q.autocalibration.transducer='transducer kunstoor - TODO';
p.calibration_profile = 'outs-2';
calibration=a3calibration(p.calibration_profile, 'stimulus_calib', ...
    {'cardgain0' 'cardgain1'}, 60, q);

p.experiment_file=[p.experiment_file_prefix '_' l num2str(p.sp_L.freq) '_' r num2str(p.sp_R.freq) '_ild' num2str(p.ild) '_' num2str(p.sp_L.modfreq) '.apx'];
stimuli=[stimuli getstim_calib()];

%%%%%%%%%%%%%%%%%%%%%%% trials %%%%%%%%%%%%%%%%%%%%%%%%%%%
trials=wraptag('trials', trials);

%%%%%%%%%%%%%%%%%%%%%%% procedure %%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.procedure
    case 2 % constant procedure
        params = [];
        params.presentations= 3;            % <presentations>3</presentations>
        params.skip         = 0;
        params.order        = 'random';     % <order>random</order>
        params.choices      = 1;            %<choices>1</choices>
        type                = p.procedure;
    otherwise
        error('Continue working on this script')
end
procedure = a3procedure(params, trials, type) ;

%%%%%%%%%%%%%%%%%%%%%%% datablocks %%%%%%%%%%%%%%%%%%%%%%%%%%%
% datablocks=[ '<datablocks>' lf '<uri_prefix>file:' xmlfilename('.') '</uri_prefix>' lf datablocks];
% datablocks=[datablocks '</datablocks>'];
datablocks=wraptag('datablocks', datablocks);


%%%%%%%%%%%%%%%%%%%%%%% stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%
stimuli=['<stimuli>' lf '<fixed_parameters>' '<parameter id="itd"/>' '</fixed_parameters>' lf stimuli lf '</stimuli>'];

%%%%%%%%%%%%%%%%%%%%%%%% corrector %%%%%%%%%%%%%%%%%%%%%%%%%%%
corrector=a3corrector('isequal'); 

%%%%%%%%%%%%%%%%%%%%%%% devices %%%%%%%%%%%%%%%%%%%%%%%%%%%

result=['<device id="soundcard" xsi:type="apex:wavDeviceType">'];
result=[result '<driver>portaudio</driver>' lf];
result=[result '<card>default</card>' lf];
result=[result '<channels>2</channels>' lf];
for count=[0:1]
    result=[result sprintf('<gain id="cardgain%d" channel="%d">0</gain>',count,count)  lf];
end
result=[result '<samplerate>44100</samplerate>' lf];
%result=[result '<blocksize>' num2str(blocksize) '</blocksize>' lf];
result=[result '<padzero>1</padzero>' lf];
result=[result '</device>' lf];
devices=wraptag('devices',result);

% devices= readfile('a3easdevices_footer.xml')];
% devices=['<devices>'];
% devices=[devices a3wavdevice('soundcard',1)];
% devices=[devices '</devices>'];
%=======================
filters=['<filters>'];
filters=[filters a3amplifier(0,'amplifier','soundcard')]; %#ok<NODEF>
filters=[filters '</filters>'];
%=======================

%%%%%%%%%%%%%%%%%%%%%%% calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%
q.autocalibration.transducer='transducer kunstoor - TODO';
calibration=a3calibration(p.calibration_profile, 'stimulus_calib', ...
    {'cardgain0' 'cardgain1'}, 60, q);

%%%%%%%%%%%%%%%%%%%%%%% connections %%%%%%%%%%%%%%%%%%%%%%%%%%%
connections=['<connections>' lf connections lf ];
connections=[connections a3connection('_ALL_',0,'soundcard', 0) lf ];
connections=[connections a3connection('datablock_calib',0,'soundcard', 1) lf ];
connections=[connections a3connection('silence',1,'soundcard',1) lf ];
connections=[connections a3connection('amplifier',0,'soundcard',1) lf];
connections=[connections '</connections>'];


%%%%%%%%%%%%%%%%%%%%%%% results %%%%%%%%%%%%%%%%%%%%%%%%%%%
results=a3results('apexresult.xsl', 'a3itd_lr_results.m','',0,1);

%%%%%%%%%%%%%%%%%%%%%%% general %%%%%%%%%%%%%%%%%%%%%%%%%%%
general=a3general(1);

%%%%%%%%%%%%%%%%%%%%%%% screens %%%%%%%%%%%%%%%%%%%%%%%%%%%
% blabels={};
% for i=1:p.intervals
%     blabels{i}=num2str(i);
% end
% screens=a3buttonscreen(blabels,'Tijdens welk interval hoor je de biep?','screen',blabels);
% screens=wraptag('screens',screens);
% if (p.intervals~=3)
%     error('Not implemented');
% end
fid=fopen([p.targetpath p.experiment_file], 'w');

il_mkexpfile;

function il_mkexpfile
    
    fwrite(fid, a3header() );
    fwrite(fid, procedure);
    fwrite(fid, corrector); %not alternatives, but isequal
    fwrite(fid, readfile('Template/ITDac_screens.xml'));
    %        fwrite(fid, screens);
    fwrite(fid, datablocks);
    fwrite(fid, devices);
    fwrite(fid, filters);
    fwrite(fid, stimuli);
    fwrite(fid, connections);
    fwrite(fid, calibration);
    fwrite(fid, results);
    fwrite(fid, general);
    fwrite(fid, a3footer());
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=itdstring(str, itd)
        if (itd==-0), itd=0; end;
    result=sprintf('%s%g', str,itd);
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=getstim_calib()
    s=['<stimulus id="stimulus_calib' '">'];
    s=[s lf '<datablocks>'];
    s=[s lf '<datablock id="datablock_calib" />'];
    s=[s lf '</datablocks>'];
    s=[s lf '<fixedParameters>'];
    s=[s lf '<parameter id="itd">100</parameter>'];
    s=[s lf '</fixedParameters>'];
    s=[s lf '</stimulus>'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_shift(insig,shift,fs)

    nsamp=round(shift*fs);
    if (abs(shift/fs)>1e6)     % microsecond accuracy
        warning(['p.shift rounded to nearest integer number of samples: ' num2str(nsamp) ]);
    end
    outsig=[zeros(nsamp,1); insig];

end

end
