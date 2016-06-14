function Create_Piano_ICRA_multi_20160428_WAE(experiment_nr,do_skip,bAdaptive)
% function Create_Piano_ICRA_multi_20160428_WAE(experiment_nr,do_skip,bAdaptive)
%
% 1. Description:
%       Generate APEX experiments if do_skip is set to 0 (default) the audio
%       signals will be generated, otherwise only the APEX (*.apx) experiments.
%       The APEX experiments are based on the template 'piano_1_2_TEMPLATE.xml'
% 
% 2. Stand-alone example:
%   Experiments = [1 1.1 2 2.1  11 11.1 12 12.1 21 21.1 22 22.1];
%   Create_Piano_ICRA_multi_20160418(Experiments);
%   
%   Experiments = [1 1.1 2 2.1]; % only C2
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
%   Experiments = 17:0.1:17.5;
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
%   Experiments = [17 17.1]; 
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 28/04/2016
% Last update on: 28/04/2016 
% Last use on   : 28/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    do_skip = 1; % 1 = keeps already created audio's
end

if nargin < 3
    bAdaptive = 1;
end

if nargin == 0
    experiment_nr = [17 17.1];
end

dir_main   = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA-v2' delim];
dir_stimuli= [dir_main 'Stimuli' delim];
dir_dst    = ['D:\Databases\dir04-Psychoacoustics\WAE\example_eval\stimuli-ABC-new' delim];
Mkdir(dir_dst);

%           0        1       2          3            4         5       6      7
pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS50','JBS73','NS19'};
 
rampout_len = 100; % ms
SNR4pede    = 20; 
note = 'C4'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parameters of the stimuli:
% bL = 1; % if 1 = then the manual set L-length is going to be used.
% fs_theo    = 44100;

dur_piano_samples = 1.5;

% L          = round(dur_piano_samples*fs_theo); % assuming 800 ms stimuli length
% loud2use   = 36; % sones
% suffixloud = sprintf('-%.0f-sone',loud2use);
% i_noises   = 12;
% gain4pede_ON = 0;
% gain4pede_OFF = -99;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
files = Get_filenames(dir_stimuli,[note 'P*.wav']);
% (1,1): C4P0t2.wav
% (1,2): C4P2t1.wav
% (1,3): C4P2t2.wav
% (1,4): C4P4t3.wav
% (1,5): C4P4t4.wav
% (1,6): C4P6t3.wav
idx_sel = [1 2 3 4 6];
files = files(idx_sel);

pairs = Get_pairwise_combinations(1,length(files));
SNR = 0;
dBnoise = 40;

pairs3= Get_triadic_combinations(1,5);

for i = 1:length(files)

    idxstr = 4;
    
    dirn = [Get_TUe_data_paths('piano') '04-PAPA\04-Background-noise\' note delim];
    
    if strcmp( files{i}(idxstr-1),'P' )
        pianoID = str2num(files{i}(idxstr));
        idx = pianoID+1;
        fname = [dirn pianos{idx} '-' note '.wav'];
        fnamenew{i} = [dirn pianos{idx} '-' note '-loop.wav'];
        [x fs] = Wavread(fname);
        
        option = [];
        option.nAnalyser = 1; % FFT
        option.CalMethod = 1;
        option.bPlot = 0;
        out = PsySoundCL(fname,option);
        Rx  = out.Data2(3:end);
        f   = out.f(3:end);
        % figure;
        % plot(f,Rx,'k-.','LineWidth',2);
        
        B = il_rx_piano(f,Rx,fs); 
        
        noise = AM_random_noise(0,fs/2,60,2,fs);
        noise = filter(B,1,noise);
        noise = noise(length(B)/2:end);
        noise = noise(1:round(dur_piano_samples*fs));
        noise = setdbspl(noise,dBnoise);
        
        % sound(noise,fs);
        Wavwrite(noise,fs,fnamenew{i});
    end  
end

for i = 1:size(pairs3,1)
    
    tmp1 = strsplit(files{pairs3(i,1)},'.');
    prefix1= tmp1{1};
    tmp2 = strsplit(files{pairs3(i,2)},'.');
    prefix2= tmp2{1};
    tmp3 = strsplit(files{pairs3(i,3)},'.');
    prefix3= tmp3{1};
    
    [insig1 fs] = Wavread([dir_stimuli files{pairs3(i,1)}]);
    [insig2   ] = Wavread([dir_stimuli files{pairs3(i,2)}]);
    [insig3   ] = Wavread([dir_stimuli files{pairs3(i,3)}]);
    
    [noise1  ] = Wavread([dir_stimuli 'rawnoises' delim 'noise_' files{pairs3(i,1)}]);
    [noise2  ] = Wavread([dir_stimuli 'rawnoises' delim 'noise_' files{pairs3(i,2)}]);
    [noise3  ] = Wavread([dir_stimuli 'rawnoises' delim 'noise_' files{pairs3(i,3)}]);
    noise      = (noise1+noise2+noise3)/3;
    
    [nback1  ] = Wavread(fnamenew{pairs3(i,1)});
    [nback2  ] = Wavread(fnamenew{pairs3(i,2)});
    [nback3  ] = Wavread(fnamenew{pairs3(i,3)});
    
    noise_backg= (nback1 + nback2 + nback3)/3;
    
    outsig1 = Do_cos_ramp(insig1,fs,0,rampout_len);
    outsig2 = Do_cos_ramp(insig2,fs,0,rampout_len);
    outsig3 = Do_cos_ramp(insig3,fs,0,rampout_len);
    
    pair_str = sprintf('%.0f%.0f%.0f',pairs3(i,:));
    Wavwrite(outsig1+noise_backg,fs,[dir_dst num2str(pairs3(i,1)) '_' prefix1 '-pair' pair_str '.wav']);
    Wavwrite(outsig2+noise_backg,fs,[dir_dst num2str(pairs3(i,2)) '_' prefix2 '-pair' pair_str '.wav']);
    Wavwrite(outsig3+noise_backg,fs,[dir_dst num2str(pairs3(i,3)) '_' prefix3 '-pair' pair_str '.wav']);
    
    Wavwrite(outsig1+noise,fs,[dir_dst num2str(pairs3(i,1)) '_' prefix1 '-pair' pair_str '-SNR-' num2str(SNR) '-dB.wav']);
    Wavwrite(outsig2+noise,fs,[dir_dst num2str(pairs3(i,2)) '_' prefix2 '-pair' pair_str '-SNR-' num2str(SNR) '-dB.wav']);
    Wavwrite(outsig3+noise,fs,[dir_dst num2str(pairs3(i,3)) '_' prefix3 '-pair' pair_str '-SNR-' num2str(SNR) '-dB.wav']);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = il_rx_piano(f,Rx,fs)

f = [0 f fs/2];
Rx = [Rx(1) Rx Rx(end)];

N = 2^12;

B = fir2(2^12,f/(fs/2),From_dB(Rx));
[H1,Fn]=freqz(B,1);

% figure;
% semilogx(Fn, 20*log10(abs([H1])));
% 
% xlabel('Frequency [Hz]')
% legend([num2str(N) ' taps']);
% title('FIR filter to be used as approximation to isolation curve')
