function [Pxx_average, F ,Pxx] = Get_PSD_analysis_arrays(insig,fs,resolution,threshold, stereoflag);
% function [Pxx_average, F ,Pxx] = Get_PSD_analysis_arrays(insig,fs,resolution,threshold, stereoflag);
%
% 1. Description:
%   Script for calculating the average Power Spectrum Density of a range of
%   wav files
%
%   Pxx_average = AnalyzePSDofWAV(filename,resolution,threshold);
%   
%   insig -
% 
%   resolution  The (minimal) frequency resolution of the spectrum in Hz
%
%   stereoflag  Determines how to cope with stereo channels (should only be
%               specified for stereo channels.
%               0 ->    left channel is analyzed
%               1 ->    the two channels are mixed before analysis
%               2 ->    the right channel is analyzed
%
% Original file name: AnalysePSDofWAV.m
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

if nargin < 3
    resolution = 20;    %The (minimal) frequency resolution of the spectrum in Hz
end
if nargin < 4
    threshold = 0.001;
end

%Derived parameters
NFFT    = 2^ceil(log(fs/resolution)/log(2)); %The number of points in the FFT
wind    = boxcar(NFFT);                      %The window used for the FFT
N       = size(insig,2);
Pxx     = zeros(NFFT/2+1,N);                 %The target matrix for all the 

for i=1:N
    %For every wav file do
    audio = insig(:,i);
    
    %Remove all silence portions out of the wav file
    if threshold > 0
        audio = removegap(audio,fs,threshold);    
    end
    
    % [Pxx(:,i),Ft] = psd(audio,NFFT,fs,wind,0,'none');%Calculate the average Power Spectrum Density of the wav file
    
    % 0 = samples for overlap, 'none' is the confidence interval
    [Pxxtmp,F] = pwelch(audio,wind,0,NFFT,fs);%Calculate the average Power Spectrum Density of the wav file
    Pxx(:,i) = Pxxtmp * (fs/2);
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