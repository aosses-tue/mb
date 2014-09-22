function r20140919_update(options)
% function r20140919_update(options)
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
% Created on    : 17/09/2014
% Last update on: 17/09/2014 % Update this date manually
% Last use on   : 17/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Voice of the Dragon

opts = Get_VoD_params;
% Voice of the dragon
ceff = 310;
L = 0.7;
n = 2:5;
fn = n*ceff/(2*L);

var2latex([n' fn' opts.mf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all

bCreateWhitenoise = 0; % Whitenoise for LIST-f

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
options.nAnalyser = 10;

bDiary = 0;
Diary(mfilename,bDiary);

tic

tmp = Get_TUe_subpaths('db_speechmaterials');
root_folder = tmp.allfiles_LISTf;
filename = [root_folder 'wdz2.wav']; % 'Elke zaterdag ga ik naar de markt'
Wavread(filename);
filenoise1 = [root_folder 'whitenoise-LISTf.wav'];
filenoise2 = [root_folder 'wivineruis.wav'];

options.CalMethod = 5;
options.bPlot = 0;

if bCreateWhitenoise
    [x1 fs1] = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
    x1RMS = rmsdb(x1);
    [x2 fs2] = Wavread([tmp.allfiles_PB 'whitenoise.wav']); % originally -10.78 dBFS RMS
    x2 = setdbspl(x2,x1RMS+100);
    rmsdb(x2)
    Wavwrite(x2,fs1,[tmp.allfiles_LISTf 'whitenoise-LISTf']);
end

out_file1 = PsySoundCL(filename,options);
out_noise1 = PsySoundCL(filenoise1,options);
out_noise2 = PsySoundCL(filenoise2,options);
f = out_file1.f;

figure;
semilogx(   f, out_file1.DataSpecOneThirdAvg, 'bo--', ...
            f, out_noise1.DataSpecOneThirdAvg-5, 'kx--', ...
            f, out_noise2.DataSpecOneThirdAvg-5,'r>-');
grid on
legend('wdz2','white noise, 60 dB SPL','SSN, 60 dB SPL')
xlabel('Frequency [Hz]')
ylabel('Sound Pressure Level [dB]')

% Saveas(gcf, 'figure','epsc') % As sent to TF, e-mail 18-09-2014

disp('');

% outputfilename = {'vod_ac_2', 'vod_ac_5'};
% lvl = 60; % dB SPL, or dB(A) if options.calmethod is set to 2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% options = Ensure_field(options,'bSave',0);
% options.calmethod   = 2; % 0 = AMT; 1 = dB(A)
% options.bPlot       = 1;
% options.modes2check = [2 5];
% options.time2save   = 5;
% options.N_periods2analyse = 3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spl2lu(20:10:90);
% h = gcf;
% 
% if options.bSave == 1
%     Saveas(h, [options.dest_folder_fig 'fig-spl-LU']);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Add window
% if options.bSave == 1 & bGenerate == 0;
%     options.bSave = ~options.bSave;
%     out = VoD_write_aligned(options); % We do not want to generate again wav-files
%     options.bSave = ~options.bSave;
% else
%     out = VoD_write_aligned(options);
% end
% 
% ti2 = options.N_periods2analyse*out.Tmodel(1);
% tf2 = ti2+(1+options.N_periods2analyse)*out.Tmodel(1);
% set( out.ha(1),'XLim',[ti2 tf2]);
% fprintf('ti = %.3f, tf = %.3f',ti2,tf2); % Use this to adapt Praat times
% out.hFig(1)
% 
% ti5 = options.N_periods2analyse*out.Tmodel(4);
% tf5 = ti5+(1+options.N_periods2analyse)*out.Tmodel(4);
% set( out.ha(2),'XLim',[ti5 tf5]);
% fprintf('ti = %.3f, tf = %.3f',ti5,tf5); % Use this to adapt Praat times
% 
% if options.bSave
%     for i = 1:length(out.hFig)
%         % options.format = 'emf'; % for PPT
%         % Saveas(out.hFig(j)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
%         options.format = 'epsc';
%         Saveas(out.hFig(i)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if bGenerate
%     
%     for i = 1:length(options.modes2check)
% 
%         mode_idx = options.modes2check(i)-1;
%         [x1 fs] = Wavread(out.Excerpt_m{mode_idx});
%         x1 = Do_calibration_level(options.calmethod,x1,lvl);
%         if options.bSave
%             disp(['Overwrite with calibrated signal (across modes), adjusted to ' Num2str(lvl) 'dBSPL'])
%             Wavwrite(x1,fs,out.Excerpt_m{mode_idx}); % Overwrite with calibrated signal
%             [bSuccess bMessage] = movefile([out.Excerpt_m{mode_idx} '.wav'], options.dest_folder);
%             if bSuccess == 0
%                 disp('File not moved, because: ')
%                 disp(bMessage)
%             end
%         end
% 
%         sil     = Gen_silence(2,fs);
% 
%         [x2 fs] = Wavread(out.Excerpt_p{mode_idx});
%         x2 = Do_calibration_level(options.calmethod,x2,lvl); %setdbspl(x2,lvl);
%         if options.bSave
%             Wavwrite(x2,fs,out.Excerpt_p{mode_idx});
%             [bSuccess bMessage] = movefile([out.Excerpt_p{mode_idx} '.wav'], options.dest_folder);
%             if bSuccess == 0
%                 disp('File not moved, because: ')
%                 disp(bMessage)
%             end
%         end
% 
%         Exp1 = ['y' num2str(i) ' = [Do_cos_ramp(x1); sil; Do_cos_ramp(x2)];'];
%         eval(Exp1);
% 
%     end
%     
%     % File to be used in Praat
%     
%     y = [y1; sil; y2];
%     if options.bSave
%         Wavwrite(y,fs,[options.dest_folder 'VoD_ac_2_and_5.wav']);
%     end
%     
% end
% 
% % fi1 = [dest_folder outputfilename{1} '.wav'];
% % fi2 = [dest_folder outputfilename{2} '.wav'];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Acoustic mode 5, modelled with-without Doppler
% t1 = ti5;
% t2 = tf5;
% fi1 = [options.dest_folder 'model-ac-mode-5-no-doppler.wav']; % predicted
% fi2 = [options.dest_folder 'model-ac-mode-5.wav']; % predicted
% options.label1 = 'Model, no doppler,';
% options.label2 = 'Model';
% [x1 fs] = Wavread(fi1);
% [x2   ] = Wavread(fi2);
% 
% if abs( rmsdb(x1) - rmsdb(x2) ) > 3
%     yxfi1 = setdbspl(x1,rmsdb(x2)+100);
%     Wavwrite(yxfi1,fs,fi2);
% end
% 
% Get_LPC_frames(x1(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-meas-5-formants-no-doppler']);
% end
% 
% options.tanalysis = [t1 t2];
% options.label = 'ac-5-doppler';
% options.LineWidth = [2 1];
% options.SPLrange = [10 70];
% 
% Get_VoD_analysis(fi1,fi2,options)
% 
% % Acoustic mode 5, measured and modelled
% fi1 = [options.dest_folder 'meas-ac-mode-5.wav']; % measured
% fi2 = [options.dest_folder 'model-ac-mode-5.wav']; % predicted
% options.label1 = ' meas';
% options.label2 = 'model';
% 
% [x1 fs] = Wavread(fi1);
% [x2 fs] = Wavread(fi2);
% 
% Get_LPC_frames(x1(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-meas-5-formants']);
% end
% 
% Get_LPC_frames(x2(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-model-5-formants']);
% end
% 
% options.tanalysis = [t1 t2];
% options.label = 'ac-5';
% Get_VoD_analysis(fi1,fi2,options);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Acoustic mode 2, modelled with-without Doppler
% t1 = ti2;
% t2 = tf2;
% fi1 = [options.dest_folder 'model-ac-mode-2-no-doppler.wav']; % predicted
% fi2 = [options.dest_folder 'model-ac-mode-2.wav']; % predicted
% options.label1 = 'Model, no doppler,';
% options.label2 = 'Model';
% [x1 fs] = Wavread(fi1);
% [x2   ] = Wavread(fi2);
% 
% if abs( rmsdb(x1) - rmsdb(x2) ) > 3
%     yxfi1 = setdbspl(x1,rmsdb(x2)+100);
%     Wavwrite(yxfi1,fs,fi2);
% end
% 
% Get_LPC_frames(x1(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-meas-2-formants-no-doppler']);
% end
% 
% options.tanalysis = [t1 t2];
% options.label = 'ac-2-doppler';
% options.LineWidth = [2 1];
% options.SPLrange = [10 70];
% 
% Get_VoD_analysis(fi1,fi2,options)% Acoustic mode 2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Acoustic mode 2
% t1 = ti2;
% t2 = tf2;
% fi1 = [options.dest_folder 'meas-ac-mode-2.wav']; % measured
% fi2 = [options.dest_folder 'model-ac-mode-2.wav']; % predicted
% options.label1 = ' meas';
% options.label2 = 'model';
% 
% [x1 fs] = Wavread(fi1);
% [x2 fs] = Wavread(fi2);
% % [flpc1 out1] = Get_LPC(x1(fs:2*fs),fs);
% % [flpc2 out2] = Get_LPC(x2(fs:2*fs),fs);
% 
% Get_LPC_frames(x1(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-meas-2-formants']);
% end
% 
% Get_LPC_frames(x2(fs:2*fs),fs);
% if options.bSave
%     h = gcf;
%     Saveas(h,[options.dest_folder_fig 'model-ac-model-2-formants']);
% end
% 
% % figure;
% % plot(   out1.f,out1.h_dB, ...
% %         out2.f,out2.h_dB )
% % xlabel('Frequency [Hz]')
% % ylabel('Gain [dB]')
% % grid on, hold on
% 
% options.tanalysis = [t1 t2];
% options.label = 'ac-2';
% options.LineWidth = [2 1];
% options.SPLrange = [10 70];
% 
% Get_VoD_analysis(fi1,fi2,options)

toc

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
