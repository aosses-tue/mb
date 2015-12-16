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
audio = wavread('C:\temp\Woman_noise-17-11-1999-CD-left channel2.wav');      %From CD woman_noise 17-11-1999
refnoise = meanrms(audio(10e5:40e5,1),44100);
clear audio;
gainSentenceNoise = ref/refnoise*10^(-5/20);    %Reference level / own level - 5dB

%Generate The Tracks for the Sentence Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indir = 'C:\Temp\SpeechTest\'; %Directory with the source material
outdir = 'C:\temp\Tracks\';
gap = 6;    %Gap in between two sentences (in seconds)
lead = 4;   %Leading noise without sentence at the beginning of each track (in seconds)
tail = 4;   %Tailing noise wihout sentence at end of track (in seconds)
fs = 44100; %Sampling frequency

%Obtain the list contents
for i=1:35  %for each list
    FID = fopen(['C:\temp\SpeechTest\Zinnen\Lijsten\lijst' num2str(i) '.opn'],'r');
    line = fgetl(FID);
    while ~strcmp(line,'path=c:\Speechtest\zinnen\wivine')
        line = fgetl(FID);
    end
    for j=1:10  %for each waveform (number/wavfile)
        line = fgetl(FID);
        ZinnenLists{i,j} = line(10:end);
    end
    fclose(FID);
end
%Lists now contains all the wav files for each list
% lists: list * item

%Generate the tracks
for i = 1:35    %for each track/list
    
    %Generate the speech channel
    s = zeros(lead*fs,1);    %add the leading silent gap
    for j = 1:10        %For each sentence
        audio = wavread([indir 'Zinnen\Wivine\' ZinnenLists{i,j}]);    %load in the sentence
        gapaudio = zeros(round(gap*fs),1);
        s = [s; audio; gapaudio];   %append the sentence and the gap to the existing track
    end
    %s = [s; zeros(round((tail-gap)*fs),1)];  %append the tailing silence to the track
    
    %Generate the noise channel (identical for each track)
    x = wavread('C:\temp\Woman_noise-17-11-1999-CD-left channel.wav',1e5+i*4e5+[1 length(s)]);  %read in the noise file
            %The noise is read in from the left channel of the noise on the CD
            %The part ead is exactly as long as the length of the speech
            %channel.
            %Every track reads in a different part from the noise track of
            %the CD
   
    %Apply the window to the noise
    x = lingate(x,0.5,0.5,fs);
    
    track = gainSentenceNoise*[s x]; %join the left and right channel
        %Apply the gain to put at -25.5 dB
    
    wavwrite(track,fs,16,[outdir 'track' num2str(i,'%.2d') '.wav']);
end
