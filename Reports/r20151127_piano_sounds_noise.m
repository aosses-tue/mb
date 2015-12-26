function r20151127_piano_sounds_noise
% function r20151127_piano_sounds_noise
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/11/2015
% Last update on: 24/11/2015 
% Last use on   : 22/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

bSave = 0;

bDoResampling   = 1; % preparing piano samples
bDoPianoNoise_example  = 0;
bDoPianoNoise   = 0;
bDoAtt4experiments = 0;
bResults = 1;

dir = [Get_TUe_data_paths('piano') '04-PAPA' delim];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoResampling
    
    notetmp = { 'Dsh', 1; ...
                'F'  , 1; ... 
                'C'  , 2; ...
                'Ash', 2; ...
                'F'  , 3;...
                'C'  , 4;... 
                'A'  , 4;...
                'Csh', 5;... 
                'C'  , 6;... 
                'G'  , 6}; 
	fstarget = 44100;
    
    for i = 1:length(notetmp)
        
        ntmp.note   = notetmp{i,1};
        ntmp.octave = notetmp{i,2};
        
        noteS    = [ntmp.note num2str(ntmp.octave)];
        f0target = note2freq(ntmp);
                
        dirtarget{i} = sprintf('%s01-Sounds%s%s%snorm-%.0f-Hz%s',dir,delim,noteS,delim,f0target,delim); % resampled-at-44100-Hz
        Resample2fs(dirtarget{i},fstarget);
        
        if i == 1
            Mkdir([dir '02-Tuned-at-44100-Hz-new'])
        end
        dstfolder = [dir '02-Tuned-at-44100-Hz-new' delim noteS delim];
        
        movefile([dirtarget{i} 'resampled-at-44100-Hz' delim],dstfolder);
        
        disp('')
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sens = 50e-3;
G    = 5;
Cal  = 1; %1/(G*sens);  % 1 =  94 dB
% Cal  = Cal/2;       % 1 = 100 dB (AMT convention)

if bDoPianoNoise_example

    dir44100 = [dir '01-Tuned-at-44100-Hz' delim];
    f1 = [dir44100 'A4' delim 'GH05-A4.wav']; 

    [x1 fs] = Wavread(f1);
    x1 = Cal*x1;
    
    Create_piano_noise(x1,fs); % get plot
    [env noise w] = Create_piano_noise(x1,fs);
    
    disp('')
end

if bDoPianoNoise
    
   dirgral  = [dir '03-Exported-as-segments' delim];
   dirgraln = [dir '04-Background-noise' delim]; Mkdir(dirgraln);
   dirs = {['C2' delim]; ['A4' delim]; ['Csh5' delim]};
   targetdur = 1; % s
   
   for i = 1:length(dirs)
       filenames = Get_filenames([dirgral dirs{i}],'.wav');
       
       for j = 1:length(filenames)
           
            [insig fs] = Wavread([dirgral dirs{i} filenames{j}]);
            Ntarget = round(targetdur*fs);
            if bSave == 1
                Wavwrite(insig(1:Ntarget),fs,[dirgraln dirs{i} num2str(targetdur*1000) '-ms-' filenames{j}]); % to be stored in noise folder
            end
            insig = Cal*insig;

            % Create_piano_noise(insig,fs); % get plot
            method = 1;
            [env noise w] = Create_piano_noise(insig,fs,method);
            
            if bSave == 1
                if j == 1
                    Mkdir([dirgraln dirs{i}]);
                end
                Wavwrite(noise,fs,[dirgraln dirs{i} 'noise-' filenames{j}]);
                Wavwrite(noise(1:Ntarget),fs,[dirgraln dirs{i} 'noise-' num2str(targetdur*1000) '-ms-' filenames{j}]);
            end        
            
       end
   end
    
end

% dir = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\02-Exported-as-segments\Csh5-mod\';
% exp4filter = 'JBS51*.wav';
% exp4filter = 'NS19*.wav';
% Get_TVL_cmd(dir,exp4filter);

if bDoAtt4experiments
    
    dirmain = ['D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot\Stimuli-A4' delim];
    dir = [dirmain 'Original' delim];
    
    fnames = Get_filenames(dir,'*ms*.wav');
    
    % % Expected list of files:
    % 1000-ms-GRAF28-A4_3.wav
    % 1000-ms-JBS51-4544-A4_4.wav
    % 1000-ms-NS19-A4_2.wav
    % noise-1000-ms-GRAF28-A4_3.wav
    % noise-1000-ms-JBS51-4544-A4_4.wav
    % noise-1000-ms-NS19-A4_2.wav
    
    Att = [0 -2 -4 0+8.64 -2 -4];
    
    for i = 1:length(Att)
        
        [x fs] = Wavread([dir fnames{i}]);
        y = From_dB(Att(i))*x;
        if bSave
            Wavwrite(y,fs,[dirmain fnames{i}]);
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dirmain = ['D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot\Stimuli-Csh5' delim];
    dir = [dirmain 'Original' delim];
    
    fnames = Get_filenames(dir,'*ms*.wav');
    
    % % Expected list of files:
    % 1000-ms-GRAF28-Cd5_3.wav
    % 1000-ms-JBS51-4544-Cd5_5.wav
    % 1000-ms-NS19-Cd5_4.wav
    % noise-1000-ms-GRAF28-Cd5_3.wav
    % noise-1000-ms-JBS51-4544-Cd5_5.wav
    % noise-1000-ms-NS19-Cd5_4.wav
    
    Att = [0 -6 -10 0 -6-9.3 -10];
    
    for i = 1:length(Att)
        
        [x fs] = Wavread([dir fnames{i}]);
        y = From_dB(Att(i))*x;
        if bSave
            Wavwrite(y,fs,[dirmain fnames{i}]);
        end
        
    end
    
    dirtmp_1 = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot\Stimuli-Csh5\';
    dirtmp = [dirtmp_1 'back\'];
    fnametmp = '1000-ms-JBS51-4544-Cd5_5.wav';
    [insig fs] = Wavread([dirtmp fnametmp]);

    Perc = -1.5; 
    outsig = Do_pitch_stretch(insig,fs,Perc,'percentage');

    fullfile_new = [dirtmp_1 fnametmp];
    if bSave
        Wavwrite(outsig,fs,fullfile_new);
    end
    
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bResults
    
    dirresults = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot\Results\';
    nRev = 6;
    
    filters = {'*C2*','*A4*','*Csh5*'};
    
    for k = 1:length(filters)
        fnames = Get_filenames(dirresults,filters{k});
        for i = 1:3 % experiment 1
            stair = a3adaptiveresults([dirresults fnames{i}]);
            [r idx] = Get_mAFC_reversals(stair.procedure1);

            try
                Th1(i,k) = median(r(end-nRev+1:end));
            catch
                Th1(i,k) = median(r);
            end
        end
        for i = 4:6 % experiment 2
            stair = a3adaptiveresults([dirresults fnames{i}]);
            [r idx] = Get_mAFC_reversals(stair.procedure1);

            try
                Th2(i-3,k) = median(r(end-nRev+1:end));
            catch
                Th2(i-3,k) = median(r);
            end
        end
    end
    
end

MinVal1 = prctile(Th1,25);
MaxVal1 = prctile(Th1,75);
MinVal2 = prctile(Th2,25);
MaxVal2 = prctile(Th2,75);
Medi1 = median(Th1);
Medi2 = median(Th2);

MiV1 = Medi1-MinVal1;
MiV2 = Medi2-MinVal2;
MaV1 = MaxVal1-Medi1;
MaV2 = MaxVal2-Medi2;

figure;
errorbar([1:3]-0.02,Medi1,MiV1,MaV1,'bo-.','LineWidth',2); hold on
errorbar([1:3]+0.02,Medi2,MiV2,MaV2,'rs-','LineWidth',2); grid on
xlim([1-0.25 3+0.25])
ylim([-15 4+4])
set(gca,'XTick',1:3)
set(gca,'XTickLAbel',{'C2','A4','Csh5'})
xlabel('Register')
ylabel('SNR, 50\%-correct scores [dB]')
legend({'JBS51-4544','GRAF28'},'Location','NorthWest')

Saveas(gcf,'results-piano-pilot')

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
