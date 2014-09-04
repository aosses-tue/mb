function r20140930_FIA_TUe(options)
% function r20140930_FIA_TUe(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 0;
%       r20140930_FIA_TUe(options);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/09/2014
% Last update on: 03/09/2014 % Update this date manually
% Last use on   : 03/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all

bGenerate = 0;

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];

bDiary = 0;
Diary(mfilename,bDiary);

% outputfilename = {'vod_ac_2', 'vod_ac_5'};
lvl = 60; % dB SPL, or dB(A) if options.calmethod is set to 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = Ensure_field(options,'bSave',0);
options.calmethod   = 2; % 0 = AMT; 1 = dB(A)
options.bPlot       = 1;
options.modes2check = [2 5];
options.time2save   = 5;

% Add window
out = VoD_write_aligned(options);
ti2 = 3*out.Tmodel(1);
tf2 = ti2+4*out.Tmodel(1);
set( out.ha(1),'XLim',[ti2 tf2]);
fprintf('ti = %.3f, tf = %.3f',ti2,tf2); % Use this to adapt Praat times
out.hFig(1)

ti5 = 3*out.Tmodel(4);
tf5 = ti5+4*out.Tmodel(4);
set( out.ha(2),'XLim',[ti5 tf5]);
fprintf('ti = %.3f, tf = %.3f',ti5,tf5); % Use this to adapt Praat times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bGenerate
    
    for i = 1:length(out.hFig)
        % options.format = 'emf'; % for PPT
        % Saveas(out.hFig(j)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
        options.format = 'epsc';
        Saveas(out.hFig(i)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
    end

    for i = 1:length(options.modes2check)

        mode_idx = options.modes2check(i)-1;
        [x1 fs] = Wavread(out.Excerpt_m{mode_idx});
        x1 = Do_calibration_level(options.calmethod,x1,lvl);
        if options.bSave
            disp(['Overwrite with calibrated signal (across modes), adjusted to ' Num2str(lvl) 'dBSPL'])
            Wavwrite(x1,fs,out.Excerpt_m{mode_idx}); % Overwrite with calibrated signal
            [bSuccess bMessage] = movefile([out.Excerpt_m{mode_idx} '.wav'], options.dest_folder);
            if bSuccess == 0
                disp('File not moved, because: ')
                disp(bMessage)
            end
        end

        sil     = Gen_silence(2,fs);

        [x2 fs] = Wavread(out.Excerpt_p{mode_idx});
        x2 = Do_calibration_level(options.calmethod,x2,lvl); %setdbspl(x2,lvl);
        if options.bSave
            Wavwrite(x2,fs,out.Excerpt_p{mode_idx});
            [bSuccess bMessage] = movefile([out.Excerpt_p{mode_idx} '.wav'], options.dest_folder);
            if bSuccess == 0
                disp('File not moved, because: ')
                disp(bMessage)
            end
        end

        Exp1 = ['y' num2str(i) ' = [Do_cos_ramp(x1); sil; Do_cos_ramp(x2)];'];
        eval(Exp1);

    end
    
    % File to be used in Praat
    
    y = [y1; sil; y2];
    if options.bSave
        Wavwrite(y,fs,[options.dest_folder 'VoD_ac_2_and_5.wav']);
    end
    
end

% fi1 = [dest_folder outputfilename{1} '.wav'];
% fi2 = [dest_folder outputfilename{2} '.wav'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic mode 2
t1 = ti2;
t2 = tf2;
fi1 = [options.dest_folder 'meas-ac-mode-2.wav']; % measured
fi2 = [options.dest_folder 'model-ac-mode-2.wav']; % predicted
options.tanalysis = [t1 t2];
options.label = 'ac-2';
options.LineWidth = [2 1];
Get_VoD_analysis(fi1,fi2,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic mode 5
t1 = ti5;
t2 = tf5;
fi1 = [options.dest_folder 'meas-ac-mode-5.wav']; % measured
fi2 = [options.dest_folder 'model-ac-mode-5.wav']; % predicted
options.tanalysis = [t1 t2];
options.label = 'ac-5';
Get_VoD_analysis(fi1,fi2,options);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
