%Generate The Tracks for the Sentence Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

indir = 'C:\Temp\SpeechTest\'; %Directory with the source material
outdir = 'C:\temp\Tracks\';
gap = 4;    %Gap in between two sentences (in seconds)
lead = 5;   %Leading noise without sentence at the beginning of each track (in seconds)
tail = 5;   %Tailing noise wihout sentence at end of track (in seconds)
fs = 44100; %Sampling frequency

%Generate the noise channel (identical for each track)
x = wavread([indir 'Zinnen\Calibration\astrid_wom.wav']);  %read in the noise file
x = x(1347:433535);         %Trim start and ending (at zero crossings)
xlong = repmat(x,20,1);     %Loop the noise track 20 times

%Generate the speech channel
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


for i = 1:35    %for each track
    s = zeros(lead*fs,1);    %add the leading silent gap
    for j = 1:10        %For each sentence
        audio = wavread([indir 'Zinnen\Wivine\' ZinnenLists{i,j}]);    %load in the sentence
        gapaudio = zeros(gap*fs,1);
        s = [s; audio; gapaudio];   %append the sentence and the gap to the existing track
    end
    s = [s; zeros((tail-gap)*fs,1)];  %append the tailing silence to the track
    
    track = [s xlong([1:length(s)],1)]; %join the left and right channel
    
    wavwrite(track,fs,16,[outdir 'track' num2str(i,'%.2d') '.wav']);
end
        
        
%Generate the tracks for the number material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate the noise channel
x = wavread([indir 'Numbers\Calibration\Numbers_calibration.wav']);  %read in the noise file
x = x(:,1); %Only retain one channel
x = x(7:557294);         %Trim start and ending (at zero crossings)
xlong = repmat(x,15,1);     %Loop the noise track 15 times

%Generate the speech channel
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
            gapaudio = zeros(gap*fs,1);
            s = [s; audio; gapaudio];   %append the sentence and the gap to the existing track
        end
        s = [s; zeros((tail-gap)*fs,1)];  %append the tailing silence to the track
    
        track = [s xlong([1:length(s)],1)]; %join the left and right channel
    
        wavwrite(track,fs,16,[outdir 'track' num2str(10*(n-1)+i+35,'%.2d') '.wav']);
    end
end
            