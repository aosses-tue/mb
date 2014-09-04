function r20140829_update(options)
% function r20140829_update(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% TO DO:
%       TIME SPECIFIC VALUES (Convert2nan_if_outofrange)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014, original name: r20140825_ENSTA_workshop.m
% Renamed on    : 28/08/2014
% Last update on: 29/08/2014 % Update this date manually
% Last use on   : 29/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all
options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
Mkdir(options.dest_folder    );
Mkdir(options.dest_folder_fig);

bGenerate = 1;
bDiary = 0;
Diary(mfilename,bDiary,options.dest_folder);

outputfilename = {'vod_ac_2', 'vod_ac_5'};
lvl = 60; % dB SPL

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

if bGenerate
    
    i = 1;
    for j = 1:2:length(out.hFig)
        options.format = 'emf'; % for PPT
        Saveas(out.hFig(j)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
        options.format = 'epsc';
        Saveas(out.hFig(j)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
        options.format = 'emf'; % for PPT
        Saveas(out.hFig(j+1),[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i)) '-JND'],options);
        options.format = 'epsc';
        Saveas(out.hFig(j+1),[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i)) '-JND'],options);
        i = i+1;
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
end    

if bGenerate % To be used in Praat
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
