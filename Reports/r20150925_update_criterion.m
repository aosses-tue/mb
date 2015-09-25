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

fs = 44100;

bPart1 = 0;
bPart2 = 1; % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart3 = 0; % Deterministic

opts.bDecisionMethod = 2;
opts.sigma      = 0.95; 
opts.audio.fs   = fs;
opts.nAnalyser  = 100; % 99.1, 100 - dau1996, my template estimation
opts.MethodIntRep = 1; % 1 - my method; 2 - using casptemplate.m

opts.Reversals4avg = 10;

if bPart1
    nFig = 3;
    [masker,insig] = exp_dau1996b(fs,nFig);
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

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     = 12;

if bPart2
    opts.DurRamps   = 150; % additional cosine ramps
    opts.bUseRamp   = 1; % additional cosine ramps
    opts.bUseRampS  = 1; % additional cosine ramps
    opts.Gain4supra =   5; % dB
    opts.audio.fs   =  fs;
    
    opts.StepdB = 2; 
    opts.StepdBmin = 0.2;
    AMTControl_cl(opts);
end

opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart3
    
    fmaskerdet{1}  = [Get_TUe_paths('outputs') 'masker-det-1.wav'];
    fmaskerdet{2}  = [Get_TUe_paths('outputs') 'masker-det-2.wav'];
    fmaskerdet{3}  = [Get_TUe_paths('outputs') 'masker-det-3.wav'];
    idxstart = [0 20 100]*1e-3*fs;
    fmaskerbuf  = [Get_TUe_paths('outputs') 'masker-buf.wav'];
    fsignals{1} = [Get_TUe_paths('outputs') 'sig1.wav'];
    fsignals{2} = [Get_TUe_paths('outputs') 'sig2.wav'];
    fsignals{3} = [Get_TUe_paths('outputs') 'sig3.wav'];
    
    try
        Wavread(fmaskerdet);
        bCreate = 0;
    catch
        bCreate = 1;
    end
    
    if bCreate
        nFig = 14;
        [masker,insig] = exp_dau1996b(fs,nFig); 

        Wavwrite(masker                 ,fs,fmaskerbuf);
        for i = 1:3
            
            maskerd = [masker(size(insig,1)+1:size(insig,1)+idxstart(i)); masker(1:size(insig,1)-(idxstart(i)))];
            Wavwrite(maskerd,fs,fmaskerdet{i});
            Wavwrite(insig(:,i),fs,fsignals{i});
        end
    end
    %%%
    
    opts.DurRamps   =  0; % additional cosine ramps
    opts.Gain4supra =  10; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 0; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.Ntimes = 1;
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Deterministic thresholds:
%     for i = 1:3
%         opts.filename1 = fmaskerdet{i};
%         opts.filename2 = fsignals{i};
%         tmp = AMTControl_cl(opts);
%         th_det(i) = tmp.Threshold;
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic thresholds:
    
    opts.Ntimes = 100;
    opts.filename1 = fmaskerbuf;
    for i = 1:3
        opts.filename2 = fsignals{i};
        tmp = AMTControl_cl(opts);
        th_sto(i) = tmp.Threshold;
    end
    
    figure;
    plot([0 20 100],th_det+75,'ro-'); hold on
    plot([0 20 100],th_sto+75,'b<--','LineWidth',2); grid on
    legend('Deterministic','Stochastic')
    
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
