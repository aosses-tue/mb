function Generate_CP810_audio_files(info)
% Script Generate_CP810_audio_files(info)
% 
% Requirements:
%   KUL_Sim_F0mod.mdl open
%   Romain Peeters' map
%   p, q structures generated. Try to use normal maps, i.e. use any map except:
%       - Maria Brughmans (21 channels enabled)
%       - Jan Leys (CSPL = 65 dB, TSPL = 35 dB)
%   Make sure SNR is set to 99 = Quiet
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bSave           = 1;
bCreateNoise    = 0;

bUsePB          = 0; % if 1 = PB; if 0 = LIST-f
bUseLIST        = ~bUsePB;
bDoSimulink     = 0;

if      bUsePB
    info.source_dir = 'E:\fda_eval\wav\';
elseif  bUseLIST
    info.source_dir = 'E:\fda_eval_LIST\wav\';
end

% % If clean-speech
% info.dest_dir = 'E:\fda_eval\wav-white\';
% p.ACC.EnableMixing  = 1;
% p.ACC.MicMixGain    = 0;
% assignin('base','p',p);

if      bUsePB
    
    info.dest_dir = 'E:\fda_eval\wav-CP810-white\';
    Speech_material = 'f';
    afiles = Get_bagshaw_names(Speech_material); % uncomment this to process male sentences
    % afiles = Get_bagshaw_names('f'); % uncomment this to process female sentences
    SNR = [99 20 10 5 0 -5 -10 -15];
    
elseif  bUseLIST
    info.dest_dir = 'E:\fda_eval_LIST\wav-CP810\';
    Speech_material = 'LISTf';
    
    directory_wav_list = 'E:\fda_eval_LIST\wav\';
    
    extra.bExtension = 0;
    afiles = Get_filenames(directory_wav_list,['wdz*.wav'],extra);
%     % afiles = {'wdz2', 'wdz4', 'wdz5', 'wdz6'};
%     afiles = {  'wdz8', 'wdz9', 'wdz12', 'wdz20', 'wdz26', 'wdz28', 'wdz30', 'wdz31', 'wdz33', 'wdz37', 'wdz40', 'wdz41', 'wdz43', ...
%                 'wdz47', 'wdz48', 'wdz51', 'wdz52', 'wdz53', 'wdz55', 'wdz58', 'wdz61', 'wdz64', 'wdz65', 'wdz68', 'wdz69', 'wdz70'};
    SNR = [99 20 10 5 0 -5];
    
%     [speechmaterial,speechpath,materiallevel,noisefile_ltass,extrainfo] = nl_speech_material(p.speechmaterial);

end

if      strcmp(Speech_material,'m')
    speechRMS = -30.69;
elseif  strcmp(Speech_material,'f')
    speechRMS = -31.3565;
elseif  strcmp(Speech_material,'LISTf')
    speechRMS = -25.2; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Opens Simulink model
open_system('KUL_Sim_F0m');

p = evalin('base','p'); 
q = evalin('base','q');
q.CFG.SimName='KUL_Sim_F0m';
assignin('base','q',q);

p.F0mod.F0modEnable = 1;
p.F0mod.F0max = 400;
assignin('base','p',p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bCreateNoise == 1
    
    try
        [noise Fsnoise] = wavread('whitenoise.wav');
    catch
        error('Generate white-noise, use output filename whitenoise.wav');
    end
    dB_source = rmsdb(noise);
    
    dB_target = speechRMS;
    
    dB2att = (dB_target - dB_source);
    amp2att = From_dB(dB2att);
    noise4speech = amp2att*noise;
    
    noisefilename = [info.source_dir 'whitenoise-' Speech_material '.wav'];
    wavwrite(noise4speech,Fsnoise, noisefilename)
    
    disp([mfilename '.m: Wav file saved as ' noisefilename])
end

if bDoSimulink
    for index_i = 71:length(afiles)

        for k = 1:length(SNR)

            q.Signal.desiredSNR = SNR(k);
            q.Signal.duration = Get_wav_duration(afiles{index_i}) + 0.5;
            assignin('base','q',q);

            set_param( q.CFG.SimName, 'StopTime', num2str(q.Signal.duration));
            set_param([q.CFG.SimName '/Normalised WAV File'], 'wavFile'   , afiles{index_i})
            if bUsePB
                set_param([q.CFG.SimName '/Normalised WAV File'], 'wavFile_Noise', ['whitenoise-' Speech_material]) 
            elseif bUseLIST
                set_param([q.CFG.SimName '/Normalised WAV File'], 'wavFile_Noise', 'wivineruis') 
            end
            set_param([q.CFG.SimName '/Normalised WAV File'], 'desiredSNR', 'q.Signal.desiredSNR')
            set_param([q.CFG.SimName '/Normalised WAV File'], 'speechRMS' , num2str(speechRMS))

    %         p.ACC.EnableMixing  = 1;
    %         p.ACC.MicMixGain    = 0;
    %         assignin('base','p',p);
            % set_param([q.CFG.SimName '/DSP Signal Path/Accesory Mixing'], 'Enable Mixing'  , num2str(noiseRMS))

            sim(q.CFG.SimName );

            while ( strcmp(get_param(q.CFG.SimName,'SimulationStatus' ),'stopped')==0 );
            end

            if ~isfield(info,'dest_dir')
                info.dest_dir = uigetdir('Select a directory to store the wav files');
                info.dest_dir = [info.dest_dir delim];
            end

            if bSave
                if q.Signal.desiredSNR == 99 % Quiet
                    wavname = ['CP810' afiles{index_i}];
                else % Conditions in noise
                    wavname = ['CP810' afiles{index_i} '-' num2str(q.Signal.desiredSNR)];
                end
                wavwrite(dir_out,p.CFG.Fs,[info.dest_dir wavname]);
                disp([mfilename '.m: Generated ' wavname '.wav'])
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            names = {'Simulation_time' 'F0'};
            alldata = dataset([], [], 'VarNames', names); 
            alldata = [alldata; dataset(F0_out.time, F0_out.signals.values, 'VarNames', get(alldata, 'VarNames'))];
            % Simulation time

            if bSave
                export(alldata, ...           % this is file name 
                            'file'          , [info.dest_dir wavname '.txt'], ...
                            'Delimiter'     , '\t', ...
                            'WriteVarNames' , true);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

% Ignore the first 2300 samples of each audio file at 15659.375 Hz
% This means the first 146.9 ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end