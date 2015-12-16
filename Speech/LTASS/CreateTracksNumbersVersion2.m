clear;clc;

%Time in between the two starts of two sentences is 6 to 7 seconds in the
%Plomp sentences (mean about 6.5 seconds)
%Silent gap in between sentences is aproximately 4.3 seconds (4 to 4.5
%seconds)

% Obtain the reference output level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The output level has to be the same as the reference signals of Wouters
%Goes Classic 1 (Vlaamse Woordenlijst)

%Calculate the meanrms of the reference signals of Wouters 1
audio = wavread('C:\Temp\Tracks of Vlaamse WoordenLijst Wouters\Track-63.wav'); %sinustoon
referencerms(1) = meanrms(audio(4e5:10e5,1),44100);
audio = wavread('C:\Temp\Tracks of Vlaamse WoordenLijst Wouters\Track-64.wav'); %1/3th octave band around 1kHz
referencerms(2) = meanrms(audio(4e5:10e5,1),44100);
audio = wavread('C:\Temp\Tracks of Vlaamse WoordenLijst Wouters\Track-69.wav'); %white noise
referencerms(3) = meanrms(audio(4e5:8e5,1),44100);
audio = wavread('C:\Temp\Tracks of Vlaamse WoordenLijst Wouters\Track-70.wav'); %pink noise
referencerms(4) = meanrms(audio(4e5:8e5,1),44100);

figure;
bar(20*log10(referencerms));title('Intensity of reference signals of ''Opname Vlaamse Woordenlijst''');
ylabel('Magnitude (dB)');
set(gca,'XTick',[1:4]);
set(gca,'XTickLabel',{'sinus';'1/3th octave band';'white noise';'pink noise'});

%We take the sinus as reference energy level:
ref = referencerms(1);

% Obtain the necessary gain/attenuation for the sentence noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The output level of the sentence-noise aught to be exactly 5 dB lower than the
%level of the reference signals

%Calculate the level of the sentence noise
audio = wavread('C:\temp\Number_noise-20-10-1999-CD-left_channel2.wav');      %From Number noise 20-10-1999
refnoise = meanrms(audio(10e5:40e5,1),44100);
clear audio;
gainNumberNoise = ref/refnoise*10^(-5/20);    %Reference level / own level - 5dB

%Generate The Tracks for the Sentence Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indir = 'C:\Temp\SpeechTest\'; %Directory with the source material
outdir = 'C:\temp\Tracks2\';
int = 5;    %Inter-number time: time between the beginning of two consecutive numbers in a list (in seconds)
lead = 4;   %Leading noise without sentence at the beginning of each track (in seconds)
tail = 4;   %Tailing noise wihout sentence at end of track (in seconds)
fs = 44100; %Sampling frequency

%Obtain the list contents
speakers = {'AG';'JW';'MD';'WD'};
for n=1:4   %for each speaker
    for i=1:10  %for each list
        FID = fopen(['C:\temp\SpeechTest\Numbers\Lijsten\listnumbers' speakers{n} num2str(i,'%.2d') '.opn'],'r');
        line = fgetl(FID);
        while ~strcmp(line,'path=c:\speechtest\Numbers')
            line = fgetl(FID);
        end
        for j=1:10  %for each waveform (number/wavfile)
            line = fgetl(FID);
            lists{n,i,j} = line(10:end);
        end
        fclose(FID);
    end
end
%Lists now contains all the wav files for each list
% lists: speaker * list * item

for n=1:4 %For each speaker
    for i=1:10  %For each list or track on the CD
        s = zeros(lead*fs,1);    %add the leading silent gap
        for j=1:10  %For each item
            audio = wavread([indir 'Numbers\' lists{n,i,j}]);    %load in the number (wavfile)
            audio = audio(:,1); %Only retain one channel
            gap = int - length(audio)/fs;
            gapaudio = zeros(round(gap*fs),1);
            s = [s; audio; gapaudio];   %append the sentence and the gap to the existing track
        end
        %s = [s; zeros((tail-gap)*fs,1)];  %append the tailing silence to the track
        
        %Generate the noise channel (identical for each track)
        x = wavread('C:\temp\Number_noise-20-10-1999-CD-left_channel.wav',1e5+i*3e5+[1 length(s)]);  %read in the noise file
            %The noise is read in from the left channel of the noise on the CD
            %The part ead is exactly as long as the length of the speech
            %channel.
            %Every track reads in a different part from the noise track of
            %the CD
   
        %Apply the window to the noise
        x = lingate(x,0.5,0.5,fs);
    
        track = gainNumberNoise*[10^(10/20)*s x]; %join the left and right channel
        %Apply the gain to put at -25.5 dB
    
        wavwrite(track,fs,16,[outdir 'track' num2str(10*(n-1)+i,'%.2d') '.wav']);
    end
end