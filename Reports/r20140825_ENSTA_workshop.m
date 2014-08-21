function r20140825_ENSTA_workshop
% function r20140825_ENSTA_workshop
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
% Created on    : 20/08/2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 20/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
dst_folder = [Get_TUe_paths('lx_Presentations') '20140825-26-BATWOMAN-ENSTA' delim 'Audios' delim];
options.dst_folder_fig = [Get_TUe_paths('lx_Presentations') '20140825-26-BATWOMAN-ENSTA' delim 'Figures-tmp' delim];
Mkdir(options.dst_folder_fig);

bGenerate = 0;
bDiary = 0;
Diary(mfilename,bDiary,dst_folder);

outputfilename = {'vod_ac_2', 'vod_ac_5'};
lvl = 70; % dB SPL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bGenerate == 1
    options.bSave = 1;
end

options.bPlot = 1;
options.modes2check = [2 5];
options.time2save = 5;
bHPF = 1;
% Add window
[xx xx out] = VoD_read_aligned(bHPF,options);
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
    for j = 1:length(out.hFig)
        options.format = 'emf'; % for PPT
        Saveas(out.hFig(j),[options.dst_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(j))],options);
        options.format = 'epsc';
        Saveas(out.hFig(j),[options.dst_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(j))],options);
    end

    for i = 1:length(options.modes2check)

        mode_idx = options.modes2check(i)-1;
        [x1 fs] = Wavread(out.Excerpt_m{mode_idx});
        x1 = setdbspl(x1,lvl);
        disp(['Overwrite with calibrated signal (across modes), adjusted to ' Num2str(lvl) 'dBSPL'])
        Wavwrite(x1,fs,out.Excerpt_m{mode_idx}); % Overwrite with calibrated signal
        movefile([out.Excerpt_m{mode_idx} '.wav'], dst_folder)

        sil     = Gen_silence(2,fs);

        [x2 fs] = Wavread(out.Excerpt_p{mode_idx});
        x2 = setdbspl(x2,lvl);
        Wavwrite(x2,fs,out.Excerpt_p{mode_idx});
        movefile([out.Excerpt_p{mode_idx} '.wav'], dst_folder)

        Exp1 = ['y' num2str(i) ' = [Do_cos_ramp(x1); sil; Do_cos_ramp(x2)];'];
        eval(Exp1);

    end
end    

% [x1 fs] = Wavread([dst_folder outputfilename{1}]);
% sil     = Gen_silence(2,fs);
% [x2 fs] = Wavread([dst_folder outputfilename{2}]);
    
if bGenerate % To be used in Praat
    
    y = [y1; sil; y2];
    Wavwrite(y,fs,[dst_folder 'VoD_ac_2_and_5.wav']);
    
end

% fi1 = [dst_folder outputfilename{1} '.wav'];
% fi2 = [dst_folder outputfilename{2} '.wav'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic mode 2
t1 = ti2;
t2 = tf2;
fi1 = [dst_folder 'modus-1_v2-2filt.wav']; % measured
fi2 = [dst_folder 'modus-1-v_2filt.wav']; % predicted
options.tanalysis = [t1 t2];
options.label = 'ac-2';
Get_VoD_analysis(fi1,fi2,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acoustic mode 5
t1 = ti5;
t2 = tf5;
fi1 = [dst_folder 'modus-4_v3-2filt.wav']; % measured
fi2 = [dst_folder 'modus-4-v_2filt.wav']; % predicted
options.tanalysis = [t1 t2];
options.label = 'ac-5';
Get_VoD_analysis(fi1,fi2,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
