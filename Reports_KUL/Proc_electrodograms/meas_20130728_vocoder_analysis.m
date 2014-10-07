function meas_20130728_vocoder_analysis(subject)
% function meas_20130728_vocoder_analysis(subject)
%
% Possible Subjects:
%   AL - Anneke Lenssen
%   AO - Alejandro Osses
%   FB - Federico Bolner
%   RZ - Rick Zweedijk
%   SJ - Sofie Jansen
%
% Programmed by Alejandro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'ConvertSequence_txt2nsb.m', ...
                    'ConvertSequence_nsb2NMT.m', ...
                    'rmsdb.m', ...
                    'Plot_sequence.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    subject = 'XX';
end

dir_roved_audio = ['/home/alejandro/Documenten/Meas/Meas/Music/Loudness_balanced_stimuli/PR_Stimuli_' subject '/'];

strat = {'ACE','F0m'};

bSave = input('Press 1 if you want to save the plots: ');

if strcmp(subject,'XX') == 0
    jmax = 2;
else
    jmax = 1;
end

for j = 1:jmax
    if strcmp(subject,'XX') == 0
        list_of_waves = {['UW_LB_' strat{j} '_131_Hz'], ...
                         ['UW_LB_' strat{j} '_139_Hz'], ...
                         ['UW_LB_' strat{j} '_147_Hz'], ...
                         ['UW_LB_' strat{j} '_156_Hz']};
    else
%         list_of_waves = {['UW_131_Hz'], ...
%                          ['UW_139_Hz'], ...
%                          ['UW_147_Hz'], ...
%                          ['UW_156_Hz']};
        list_of_waves = {['UW_104_Hz'], ...
                         ['UW_110_Hz'], ...
                         ['UW_117_Hz'], ...
                         ['UW_124_Hz'], ...
                         ['UW_131_Hz']};
    end
                 
    for i = 1:length(list_of_waves)
        filename = [dir_roved_audio list_of_waves{i}];
        [x, Fs] = wavread(filename);

        deltaValue(i) = quick_check_RMS_wav(filename);
        N = 4096*4;
        [H_dB(:,i), f] = semilogx_fft(x,1,Fs,length(x));
        LegendString{i} = [name2figname(list_of_waves{i}) ' \Delta dB = ' num2str(deltaValue(i))];
    end

    figure
    semilogx(f, H_dB), grid on
    xlim([100 1000])
    legend(LegendString, 'Location', 'NorthWestOutside')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(subject,'XX') == 0
        list_of_waves = {['UW_LB_' strat{j} '_208_Hz'], ...
                     ['UW_LB_' strat{j} '_220_Hz'], ...
                     ['UW_LB_' strat{j} '_233_Hz'], ...
                     ['UW_LB_' strat{j} '_247_Hz']};
    else
        list_of_waves = {   ['UW_208_Hz'], ...
                            ['UW_220_Hz'], ...
                            ['UW_233_Hz'], ...
                            ['UW_247_Hz']};
    end

    for i = 1:length(list_of_waves)
        filename = [dir_roved_audio list_of_waves{i}];
        [x, Fs] = wavread(filename);

        deltaValue(i) = quick_check_RMS_wav(filename);
        
        [H_dB(:,i), f] = semilogx_fft(x,1,Fs);
        LegendString{i} = [name2figname(list_of_waves{i}) ' \Delta dB = ' num2str(deltaValue(i))];
    end

    figure
    semilogx(f, H_dB), grid on
    xlim([100 1000])
    legend(LegendString, 'Location', 'NorthWestOutside')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirAudioFiles = '/home/alejandro/Documenten/Meas/Meas/Music/Loudness_balanced_stimuli/Vocoder-out-without-roving/';
addpath(dirAudioFiles);

sequences_from_diet = { '-104', ...
                        '-131', ...
                        '-165', ...
                        '-208', ...
                        '-262'};
h = [];

Tmin_CP = 0:1:4;
Tmax_CP = Tmin_CP + 1;
EnergyACE = [];
EnergyF0m = [];
        
for i = 1:length(sequences_from_diet)
    for j = 1:2
        sequence_name = [strat{j} sequences_from_diet{i} '.out'];
        seq = ConvertSequence_txt2nsb(sequence_name);
        p = ConvertSequence_nsb2NMT(seq);
        Plot_sequence(p);
        title(sequence_name);
        h(end+1) = gca;
        
        Tmin_CP_idx  = getIndexElectrodogram(Tmin_CP(j), p );
        Tmax_CP_idx  = getIndexElectrodogram(Tmax_CP(j), p   );
        
        if j == 1
            EnergyACE   = [EnergyACE Average_per_Channel_per_segment( p    , Tmin_CP_idx   , Tmax_CP_idx    )]; % [CU/s]
        else
            EnergyF0m   = [EnergyF0m Average_per_Channel_per_segment( p    , Tmin_CP_idx   , Tmax_CP_idx    )];
        end
        
        wavfile_name = [strat{j} sequences_from_diet{i} '.wav'];
        quick_Spectrogram(wavfile_name);
        
        if bSave == 1
            plotfilename = [dirAudioFiles strat{j} sequences_from_diet{i}];
            saveas(gcf, plotfilename,'epsc');
            disp(['Plot saved as: ' plotfilename])
        end
        
    end
end

linkaxes(h,'xy')
xlim([0 5*1000])

figure
bar(EnergyACE), legend({'0', '1', '2', '3', '4'})
title('ACE')

figure
bar(EnergyF0m), legend({'0', '1', '2', '3', '4'})
title('F0mod')

rmpath(dirAudioFiles);

end
