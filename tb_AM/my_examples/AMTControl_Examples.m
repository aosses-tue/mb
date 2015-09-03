function filename = AMTControl_Examples(nExample)
% function filename = AMTControl_Examples(nExample)
%
% 1. Description:
%       Wav files are generated if no output is specified
%   
% 2. Stand-alone example:
%       % Tone generation + save to file of audios of example 1
%       AMTControl_Examples(1); 
%       
%       AMTControl_Examples(3); 
% 
%       % To get the filenames but not storing the audio files
%       filename = AMTControl_Examples(1); 
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20150522_update (for example 1)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 10/08/2015
% Last update on: 10/08/2015 
% Last use on   : 19/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dirout = [Get_TUe_paths('outputs') 'AMTControl-examples' delim];
dBFS = 100;

if nargin == 0
    nExample = 1;
end

switch nExample
    case 0 % default
        
        file1 = ['dau1996b_expI_noisemasker.wav'];
        file2 = ['dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav'];
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
    case 1 % Weber experiment
        
        fs = 44100;
        f = 1000;
        dur = 4; % in seconds
        SPLs = [60 42]; % see r20150522_update. Estimated dprime should be around 1.25
        
        [outsig1 file1] = Il_create_tone(f,dur,fs,SPLs(1));
        [outsig2 file2] = Il_create_tone(f,dur,fs,SPLs(2));
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
        if nargout == 0
            Wavwrite( outsig1,fs,filename{1} );
            Wavwrite( outsig2,fs,filename{2} );
        end
        
    case 2 % Fluctuation strength samples
        
        file1 = ['fluct_test_bbn_AM_m_000_fmod_004Hz_60_dBSPL.wav'];
        file2 = ['fluct_test_bbn_AM_m_070_fmod_004Hz_60_dBSPL.wav'];
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
    case 3 % Modulation detection
        
        % - An independent sample of noise was presented at each interval
        % - Noise stimuli were digitally filtered before modulation = first filtered, then modulated
        % - To eliminate level cues , digital waveforms were adjusted to have equal power
        % - Modulation with zero onset  phase was applied
        % - Both were windowed with 200-ms cosine ramps
        
        fs = 44100;
        BW = 3;
        fc = 5000;
        fmod = 20;
        dur = 10; % buffer longer than what we need...
        Mdept = 0.5; % well above threshold
        SPL = 70;
        fs = 44100;
        [outsig1 file1,  xx, outsigBBN] = AM_random_noise_BW(fc,BW,SPL,dur,fs,fmod,0);
        [xx      file2, env           ] = AM_random_noise_BW(fc,BW,SPL,dur,fs,fmod,Mdept);

        outsig2 = outsig1.*env;
        outsig2 = setdbspl(outsig2,SPL,dBFS);
        
        filename{1} = [dirout file1{1} '.wav'];
        filename{2} = [dirout file2{1} '.wav'];
        filename{3} = [dirout file1{2} '.wav'];
        
        if nargout == 0
            Wavwrite( outsig1,fs,filename{1} );
            Wavwrite( outsig2,fs,filename{2} );
            % Wavwrite( outsigBBN,fs,filename{3} );
        end
    case 4
        % 4 - Envelope detection in two auditory channels
        % Hartmann2005, 'Useless envelope', Ex 3:
        f = [100 1000];
        A = [1 1/4];
        phi = 0;
        dur = 500e-3; 
        fs = 44100;
        outsig1 = zeros(dur*fs,1);
        for i = 1:2
            ytmp = A(i) * Create_sin_phase(f(i),phi,dur,fs);
            outsig1 = outsig1+ytmp;
        end

        t = (1:length(outsig1)) /fs;
        env = sqrt( 17/16 + 0.5*cos(2*pi*900*t) );

        % figure;
        % plot(t*1000, env,'r',t*1000,outsig1,'b'); hold on
        % plot(t*1000,-env,'r')
        % xlabel('Time [ms]')
        % ylabel('Amplitude')
        % grid on
        % legend('Envelope')

        SPL = 70;
        outsig1 = setdbspl(outsig1,SPL,dBFS);
        outsig1 = [Gen_silence(50e-3,fs); outsig1];
        fname = sprintf('%ssine-%.0f-plus-%.0f-Hz-%.0f-dB.wav',dirout,f(1),f(2),SPL);
        if nargout == 0
            Wavwrite(outsig1,fs,fname);
        end
        
        filename{1} = fname;
        filename{2} = fname;
        
    case {5, 6}
        
        fs = 44100;
        durM = 10; % 500e-3;
        
        BW = 100;
        fc = 1300;
        fmod = 0;
        SPL = 60;
        
        % BP = AM_random_noise(500,800,SPL-25,durM,fs); % background noise to avoid undesired cues
        warning('BP bypassed');
        
        [outsig1 file1,  xx, outsigBBN] = AM_random_noise_BW(fc,BW,SPL,durM,fs,fmod,0);
        % outsig1 = outsig1 + BP;
        
        filename1 = [dirout file1{1} '.wav'];
        
        % Test tone:
        durT = 400e-3; % plus 50 and 50 ms of silence
        dursilence = 50e-3;
        f = 2000;
        rampupdn = 20;
        % test tone:
        outsig2 = Il_create_tone(f,durT,fs,SPL,rampupdn,dursilence);
        
        fname = sprintf('%ssine-%.0f-Hz-ramps-of-%.0f-ms-%.0f-dB.wav',dirout,f,rampupdn,SPL);
        filename2 = fname;
        
        % Multiplication noise:
        outsig3 = Multiplied_noise(fc,BW,SPL,durM,fs);
        outsig3 = setdbspl(outsig3,SPL);
        % outsig3 = outsig3 + BP;
              
        fname = sprintf('%smult-noise-%.0f-Hz-BW-%.0f-Hz-%.0f-dB.wav',dirout,fc,BW,SPL);
        filename3 = fname;
                
        if nExample == 5;
            filename{1} = filename1;
        elseif nExample == 6
            filename{1} = filename3;
        end
        filename{2} = filename2;
        
        if nargout == 0
            Wavwrite(outsig1,fs,filename1);
            Wavwrite(outsig2,fs,filename2);
            Wavwrite(outsig3,fs,filename3);
        end    
        
    case {7, 8}
        
        fs = 44100;
        durM = 10; 
        
        BW = 20;
        fc = 1300;
        fmod = 0;
        SPL = 60;
        
        %BP = AM_random_noise(500,800,SPL-25,durM,fs); % background noise to avoid undesired cues
        warning('BP bypassed');
        
        [outsig1 file1,  xx, outsigBBN] = AM_random_noise_BW(fc,BW,SPL,durM,fs,fmod,0);
        % outsig1 = outsig1 + BP;
        
        filename1 = [dirout file1{1} '.wav'];
        
        % Test tone:
        durT = 400e-3; % plus 50 and 50 ms of silence
        dursilence = 50e-3;
        f = 2000;
        rampupdn = 20;
        % test tone:
        outsig2 = Il_create_tone(f,durT,fs,SPL,rampupdn,dursilence);
        
        fname = sprintf('%ssine-%.0f-Hz-ramps-of-%.0f-ms-%.0f-dB.wav',dirout,f,rampupdn,SPL);
        filename2 = fname;
        
        % Multiplication noise:
        outsig3 = Multiplied_noise(fc,BW,SPL,durM,fs);
        outsig3 = setdbspl(outsig3,SPL);
        % outsig3 = outsig3 + BP;
              
        fname = sprintf('%smult-noise-%.0f-Hz-BW-%.0f-Hz-%.0f-dB.wav',dirout,fc,BW,SPL);
        filename3 = fname;
                
        if nExample == 7;
            filename{1} = filename1;
        elseif nExample == 8
            filename{1} = filename3;
        end
        filename{2} = filename2;
        
        if nargout == 0
            Wavwrite(outsig1,fs,filename1);
            Wavwrite(outsig2,fs,filename2);
            Wavwrite(outsig3,fs,filename3);
        end 
        
    case 9
        
        fs = 44100;
        durM = 10; % 500e-3;
        fc = 1300;
        SPL = 60;
        
        % BP = AM_random_noise(500,800,SPL-25,durM,fs); % background noise to avoid undesired cues
        warning('BP bypassed');
        
        % Sine-masker:
        outsig1 = .5*Create_sin(fc,durM,fs,0);
        outsig1 = setdbspl(outsig1,SPL);
        % outsig1 = outsig1 + BP;
        
        fname = sprintf('%ssine-%.0f-Hz-%.0f-dB.wav',dirout,fc,SPL);
        filename1 = fname;
        filename{1} = filename1;
        
        % Test tone:
        durT = 400e-3; % plus 50 and 50 ms of silence
        dursilence = 50e-3;
        f = 2000;
        rampupdn = 20;
        % test tone:
        outsig2 = Il_create_tone(f,durT,fs,SPL,rampupdn,dursilence);
        
        fname = sprintf('%ssine-%.0f-Hz-ramps-of-%.0f-ms-%.0f-dB.wav',dirout,f,rampupdn,SPL);
        filename2 = fname;
        
        filename{2} = filename2;
        
        if nargout == 0
            Wavwrite(outsig1,fs,filename1);
            Wavwrite(outsig2,fs,filename2);
        end    
        
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [outsig filename] = Il_create_tone(f,dur,fs,lvl,rampupdn,dursilence)

% Creating test tone:
y       = .5*Create_sin(f,dur,fs,0);
if nargin < 5
    rampupdn  = 5; % ms
end
if nargin < 6
    dursilence = 200e-3;
end

insig    = Do_cos_ramp(y,fs,rampupdn,rampupdn);

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(dursilence,fs); ...
               outsigtmp;
               Gen_silence(dursilence,fs)];
           
filename = sprintf('tone-f-%.0f-Hz-at-%.0f-dB-dur-%.0f-s.wav',f,lvl,dur);
