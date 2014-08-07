% function y = VoD_comparisons_do_Zwicker(x)
% function y = VoD_comparisons_do_Zwicker(x)
%
% 1. Description:
%       Only use it called from VoD_comparisons2
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:   
%       VoD_comparisons2(0);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
info.bNewFigure = 0;
info.color = 'b';

info.filename = misc.near_field_filename{mode_idx,1};

info.bSave = 0;
res{1} = Zwicker_dynamic_loudness_model(ymeas, fs, info);

info.color = 'r';
info.filename = misc.near_field_filename{mode_idx,2};
res{2} = Zwicker_dynamic_loudness_model(ymodel, fs, info);
info.bSave = 1;

subplot(2,1,1)
legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
halocal = gca;

scale = 1.8;
Ylimit = get(halocal(end),'YLim');
set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

subplot(2,1,2)
legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
halocal(end+1) = gca;
scale = 1.8;
Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

linkaxes(halocal,'x')

info = rmfield(info,'filename');
% info.bSave = 0;

h_Loudness(end+1) = gcf;

if mode_idx == 4

    ii = 1;
    figure;
    info.bNewFigure = 0;
    info.color = 'b';

    Nfinal = round(misc.Tmodel(ii)*2*44100); %  seconds

    filename1 = 'modus-1_v2-2filt-fc-251-Hz.wav';
    [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
    ymeas2 = ymeas2(1:Nfinal);
    info.bSave = 0;
    res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

    filename2 = 'modus-1-v_2filt-fc-251-Hz.wav';
    [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
    ymodel2 = ymodel2(1:Nfinal);
    info.color = 'r';
    res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
    info.bSave = 1;

    subplot(2,1,1)
    legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
    halocal = gca;
    title(name2figname(filename1))

    scale = 1.8;
    Ylimit = get(halocal(end),'YLim');
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

    subplot(2,1,2)
    legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
    halocal(end+1) = gca;
    Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);
    hZwicker_band = gcf;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High frequency band

    figure;

    filename1 = 'modus-1_v2-2filt-fc-1000-Hz.wav';
    [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
    ymeas2 = ymeas2(1:Nfinal);
    info.bSave = 0;
    info.color = 'b';
    res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

    filename2 = 'modus-1-v_2filt-fc-1000-Hz.wav';
    [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
    ymodel2 = ymodel2(1:Nfinal);
    info.color = 'r';
    res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
    info.bSave = 1;

    subplot(2,1,1)
    legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
    halocal(end+1) = gca;
    title(name2figname(filename1))

    scale = 1.8;
    Ylimit = get(halocal(end),'YLim');
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

    subplot(2,1,2)
    legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
    halocal(end+1) = gca;
    Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

    hZwicker_band(end+1) = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High frequency band

    figure;

    filename1 = 'modus-1_v2-2filt-fc-3981-Hz.wav';
    [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
    ymeas2 = ymeas2(1:Nfinal);
    info.bSave = 0;
    info.color = 'b';
    res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

    filename2 = 'modus-1-v_2filt-fc-3981-Hz.wav';
    [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
    ymodel2 = ymodel2(1:Nfinal);
    info.color = 'r';
    res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
    info.bSave = 1;

    subplot(2,1,1)
    legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
    halocal(end+1) = gca;
    title(name2figname(filename1))

    scale = 1.8;
    Ylimit = get(halocal(end),'YLim');
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

    subplot(2,1,2)
    legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
    halocal(end+1) = gca;
    Ylimit = get(halocal(end),'YLim'); 
    set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

    linkaxes(halocal,'x')

    hZwicker_band(end+1) = gcf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
