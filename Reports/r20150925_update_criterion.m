function r20150925_update_criterion
% function r20150925_update_criterion
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/09/2015
% Last update on: 22/09/2015 
% Last use on   : 22/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

fs = 22050;

bPart1 = 0;
bPart2 = 1;

opts.bDecisionMethod = 4;
opts.sigma      = 1.85; % 0.5
opts.audio.fs   = fs;
opts.nAnalyser  = 100; % 99.1, 100 - dau1996, my template estimation
opts.MethodIntRep = 1; % 1 - my method; 2 - using casptemplate.m

if bPart1
    [masker,insig] = exp_dau1996b(fs);
    fmasker  = [Get_TUe_paths('outputs') 'masker.wav'];
    fsignals{1} = [Get_TUe_paths('outputs') 'sig1.wav'];
    fsignals{2} = [Get_TUe_paths('outputs') 'sig2.wav'];
    fsignals{3} = [Get_TUe_paths('outputs') 'sig3.wav'];

    Wavwrite(masker,fs,fmasker);
    for i = 1:3
        Wavwrite(insig(:,i),fs,fsignals{i});
    end
    opts.DurRamps = 0; % additional cosine ramps
    opts.Gain4supra = 0; % dB

    opts.filename1 = fmasker;
    for i = 3
        opts.filename2 = fsignals{i};
        opts.do_template = 1;
        opts.do_simulation = 0;

        % opts.DurRamps = 150; % additional cosine ramps
        % opts.Gain4supra = 5; % dB
        
        AMTControl_cl(opts);
    end
end

if bPart2
    opts.do_template = 1;
    opts.do_simulation = 1;
    opts.DurRamps = 150; % additional cosine ramps
    opts.Gain4supra = 5; % dB
    opts.audio.fs   = fs;
    opts.Nreversals = 12;
    
    AMTControl_cl(opts);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
