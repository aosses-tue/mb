function [Pxx_average, F ,Pxx] = AnalysePSDofWAV(filename,resolution,threshold, stereoflag);
% function [Pxx_average, F ,Pxx] = AnalysePSDofWAV(filename,resolution,threshold, stereoflag);
%
% 1. Description:
%   Script for calculating the average Power Spectrum Density of a range of
%   wav files
%
%   Pxx_average = AnalyzePSDofWAV(filename,resolution,threshold);
%   
%   filename    is the filename of the wav file or the direcotry of the wav
%               files if more than 1 wav is processed
%   resolution  The (minimal) frequency resolution of the spectrum in Hz
%
%   stereoflag  Determines how to cope with stereo channels (should only be
%               specified for stereo channels.
%               0 ->    left channel is analyzed
%               1 ->    the two channels are mixed before analysis
%               2 ->    the right channel is analyzed
%
% % Example:
%       % All wav files inside the Matrix directory:
%       filename = 'D:\Databases\dir03-Speech\spanish\Matrix\';
% 
%       % Only one wav file:
%       filename = 'D:\Databases\dir03-Speech\spanish\Matrix\00131.wav';
% 
%       [filename fileplusdir] = Get_filenames('D:\Databases\dir03-Speech\spanish\Matrix\',['*.wav']);
%       filename = fileplusdir(1:260);
%       [Pxx_average, F ,Pxx] = AnalysePSDofWAV(filename);
%       figure; plot(F,Pxx_average); grid on
%       save([Get_TUe_paths('outputs') 'EsMatrixSpectra'],'F','Pxx','Pxx_average');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isunix)
    delim='/';
else
    delim='\';
end

%Check the filename(s) of the wav files
if iscell(filename) 
    %The filename is a cell array of wav files
    filenames = filename;
elseif strcmp(filename(end),delim)
    %The filename is a directory -> Read all wav files in the directory
    dirname = filename;
    d = dir([dirname '*.wav']);
    N = length(d);
    for i=1:N
        filenames{i} = [dirname d(i).name];
    end
elseif strcmp(filename(end-3:end),'.wav')
    %The filename is a wav file -> Read in just 1 filename
    filenames{1} = filename;
else
    error('The file is not supported (only wav files or entire directories).');
end

if nargin < 2
    resolution = 20;    %The (minimal) frequency resolution of the spectrum in Hz
end
if nargin < 3
    threshold = 0.001;
end

[audio,rate1] = wavread(filenames{1});
    
%Derived parameters
NFFT    = 2^ceil(log(rate1/resolution)/log(2)); %The number of points in the FFT
wind    = boxcar(NFFT);                         %The window used for the FFT
N       = length(filenames);
Pxx     = zeros(NFFT/2+1,N);                    %The target matrix for all the 

for i=1:N
    %For every wav file do
    disp(filenames{i});
    [audio,rate] = wavread(filenames{i});
    if (size(audio,2) ~= 1)
        if nargin < 4
            error('Only mono wav files supported (or select a stereoflag)');
        elseif stereoflag == 0
            audio = audio(:,1);
        elseif stereoflag == 1
            audio = mean(audio,2);
        elseif stereoflag == 2
            audio = audio(:,2);
        else
            error('Un understood stereoflag');
        end
    end
    if (rate ~= rate1)
        error('All the wav files have to have the same sampling rate');
    end
    
    %Remove all silence portions out of the wav file
    if threshold > 0
        audio = removegap(audio,rate,threshold);    
    end
    
    [Pxx(:,i),F] = psd(audio,NFFT,rate,wind,0,'none');  %Calculate the average Power Spectrum Density of the wav file
    Nsections(i) = ceil(length(audio)/NFFT);        %The number of blocks (of NFFT samples) that are in the wavfile
    Pxx(:,i) = Nsections(i)*Pxx(:,i);               %Weigth the importance of this wav files Pxx by Nsections
    
end

%Average over all wav files
Pxx_average = sum(Pxx,2)/sum(Nsections);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Helper function to remove the silenced gaps in the between words
function y = removegap(x,fs,thr);
blocklength = round(0.02*fs);  %Use a 20ms window
VA = zeros(size(x));
for i = 1:length(x)
    VA(i) = (sqrt(mean(x(i:min(i+blocklength-1,length(x))).^2)) > thr);
end
y = x(find(VA));
return;