function res = Zwicker_dynamic_loudness_model(sig,fs,info)
% function res = Zwicker_dynamic_loudness_model(sig,fs,info)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%   Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on : 30/06/2014
% Last update: 03/07/2014 % Update this date manually
% Last used  : 03/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fs < 44100
    sig = resample(sig,44100,fs);
    disp([mfilename '.m: signal resampled to 44.1 kHz, to avoid ''filter design'' (CB) problems...'])
    fch = fs/2 - 1;
    fs  = 44100;
    
    if info.bSave
    
        outputdir = [Get_TUe_paths('outputs') 'new' delim];
        % Mkdir(outputdir);
        info.fs = fs;
        sig     = freqfftwlpf(sig,info,fch);
        Wavwrite(sig,fs,[outputdir Split_file_name(info.filename)]);
        
    end
    
end

if nargin < 3
    info.bNewFigure = 1;
end

info = Ensure_field(info,'color','b');

res = Loudness_TimeVaryingSound_Zwicker(sig, fs);
res.sig = sig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[res.InstantNs_MA_1_Bark, res.InstantNs_MA] = Reduce_Bark_resolution(res);

% figure; 
% plot(res.time, sum(res.InstantNs_MA'       )*res.BarkStep,'r'), hold on
% plot(res.time, sum(res.InstantNs_MA_1_Bark')             ,'b--')

res2print = 100*sum(res.InstantNs_MA_1_Bark)/sum(sum(res.InstantNs_MA_1_Bark));
disp('Energy [%] per Critical Band')
for i = 1:length(res2print)
    fprintf('%.0f:\t %.2f [%%]\n',i,res2print(i))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if info.bNewFigure == 1
    figure;
else
    hold on;
end

t       = res.time;
N_inst  = res.InstantaneousLoudness;
tsig    = 0:1/fs:(length(sig)-1)/fs;

subplot(2,1,1);
plot(t, N_inst,info.color); hold on, grid on

res.txt_loudness = sprintf('N_t = %.2f [sone], %.1f [phon]',res.Nt,res.Lt);
legend(res.txt_loudness)

xlabel('time [s]');
ylabel('Loudness [sones]');

subplot(2,1,2)
plot(tsig, sig,info.color); hold on, grid on

res.txt_db = sprintf('Lrms = %.1f [dB]',rmsdb(sig));
legend(res.txt_db)

xlabel('time [s]');
ylabel('Pressure, Pa');

res.plot_time       = repmat(res.time'   ,                1,length(res.barkAxis));
res.plot_barkAxis   = repmat(res.barkAxis, length(res.time),                  1 ); 

if nargin == 0
    figure;
    mesh(res.plot_time, res.plot_barkAxis, res.InstantaneousSpecificLoudness)
    xlabel('Time [s]')
    ylabel('Critical-band rate [bark]')
    zlabel('Specific loudness [sone]')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end