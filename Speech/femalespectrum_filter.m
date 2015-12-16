function [B]=femalespectrum_filter(fs, type)
% function [B]=femalespectrum_filter(fs, type)
%
% 1. Description:
%
% 2. Stand-alone example:
%       femalespectrum_filter;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<1)
    fs=44100;
end
if (nargin<2)
    type='ansi';
end

cross=[200 500];
error('This is still malespectrum_filter: in adaptation process')

%F=[0:100:fs/2];                 % frequencies
F=logspace(0, log10(fs/2));
F(1)=0;
F(end)=fs/2;

switch type
    case 'ansi'
        
        R=zeros(1,length(F));             % response
        
        i=1;
        ref=log2(cross(1));
        while (F(i)<cross(1))
            R(i)=10^((ref-log2(F(i)))*(-22)/20); % -12
            i=i+1;
        end
        while (F(i)<cross(2))
            R(i)=1;
            i=i+1;
        end
        ref=log2(cross(2));
        while (i<=length(F))
            R(i)=10^((log2(F(i))-ref)*(-9)/20);
            i=i+1;
        end
        
        
    case 'moore'
%         PSD_SPAN = 1000*round(fs/1000); %% need reasonable resolution of psd, does not have to be power of 2
%         f = ((1:PSD_SPAN)-0.5)*(fs/2)/PSD_SPAN;
        bindB   = 0*F - 100; %% default to really small number
        
        %% flat portion 100-500 Hz
        span = find((F>99) & (F<=500));
        bindB(span) = 0; %% power per bin, 0 dB to start with
        corner = max(span);
        %% sloping -7.5dB/oct above 500 Hz but only up to 8 kHz : WAS -9dB/oct above 500 Hz to Nyquist
        logf = log10(F);
        near8k = find( F<=8000 , 1, 'last' );
        midspan = corner+1:near8k;
        bindB(midspan) = -7.5*( logf(midspan) - logf(corner) )/log10(2); %% implicit +0 since flat response was 0 dB
        %% NEW sloping -13dB/oct above 8 kHz
        endspan = near8k+1:length(bindB);
        bindB(endspan) = -13 *( logf(endspan) - logf(near8k) )/log10(2) + bindB(near8k); %% absolute change in response from 8k response
        
        R=10.^(bindB/20);
        
    otherwise
        error('Invalid type');
end

B=fir2(10000,F/fs*2, R);

doplot=1;
if (doplot)
    figure
    plot(F,20*log10(R),'-xr');
    title('Male filter');
    grid on;
    hold on
    [resp,Fresp]=freqz(B,1,[],fs);
    plot(Fresp,20*log10(abs(resp)), '-ob');
    
    %% Reference: SII
    f = [160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
        2500 3150 4000 5000 6300 8000];
    E = SpeechSptr('normal');
    plot(f,E-max(E),'g');
    
    legend('Target', 'Filter', 'ANSI S3.5 normal');
end

disp('')
