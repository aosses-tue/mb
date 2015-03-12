function r=vocoder(d,fs,nbands)
% function r=vocoder(d,fs,nbands)
%
% Copyright Tom Francart, 2009
% Edited by: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
    fs=44100;
end
if (nargin<3)
    nbands=8;
end

type='freedom';
doplot=0;

%% Determine analysis crossover frequencies
analysis_cutoff=analysis_cutoff_freqs_cochlear(nbands, type);

%% Create analysis filter bank
filter_order=4;
analysis_filters_B=zeros(filter_order*2+1,nbands);
analysis_filters_A=analysis_filters_B;
for i=1:nbands
    [analysis_filters_B(:,i), analysis_filters_A(:,i)] = ...
        butter(filter_order, ...
            [analysis_cutoff(i) analysis_cutoff(i+1)]/fs*2);
end

if (doplot)
    plot_filterbank(analysis_filters_B, analysis_filters_A);
    title('Analysis filters');
end

%% Create synthesis filter bank
resynthesis_filters_B=analysis_filters_B;
resynthesis_filters_A=analysis_filters_A;

% resynthesis_cutoff=resynthesis_cutoff_freqs(nbands, 'greenwood');
% for i=1:nbands
%     [resynthesis_filters_B(:,i), resynthesis_filters_A(:,i)] = ...
%         butter(filter_order, ...
%             [resynthesis_cutoff(i) resynthesis_cutoff(i+1)]/fs*2);
% end

if (doplot)
    plot_filterbank(resynthesis_filters_B, resynthesis_filters_A);
    title('Resynthesis filters');
end

%% Create LP filter for envelope detection
[Blp,Alp]=butter(4, 200/fs*2);

%% Process input signal
out=zeros(length(d),nbands);
noise=randn(length(d),1);
for i=1:nbands
    % Filter input signal
    t=filter(analysis_filters_B(:,i), analysis_filters_A(:,i), d);
    
    % Envelope detection
    t=(t+abs(t))/2;
    t=filter(Blp,Alp,t);
    
    % Create noise band
    noiseband=filter(resynthesis_filters_B(:,i), ...
        resynthesis_filters_A(:,i), noise);
    
    out(:,i)=t.*noiseband;
end

r=sum(out,2);



function plot_filterbank(B,A)

figure;
hold all;

for i=1:size(B,2)
    [H,F]=freqz( B(:,i), A(:,i), [], 44100);
    semilogx(F, 20*log10(abs(H)));
end

xlabel('Frequency');
ylabel('Magnitude');
set(gca,'XScale', 'log');


