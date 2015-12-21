function [B]=femalespectrum_filter(fs, type,bPlot)
% function [B]=femalespectrum_filter(fs, type,bPlot)
%
% 1. Description:
%       It generates a female spectrum filter considering a normal voice effort.
%       Input parameters:
%           - fs: sampling frequency [Hz]
%           - type = 'ansi'  for using band-pass filter between 100-500 Hz 
%                            with lower and upper slopes of 22 and -9 dB/oct
%       The 'ansi' version was used in Dreschler2001.
% 
% 2. Stand-alone example:
%       femalespectrum_filter;
% 
% 3. Additional information:
%       Tested cross-platform: Yes
%       See also SpeechSptr.m; malespectrum_filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<1)
    fs=44100;
end
if (nargin<2)
    type='ansi';
end
if nargin < 3
    bPlot = 0;
end

cross=[200 500];

F=logspace(0, log10(fs/2));
F(1)=0;
F(end)=fs/2;

switch type
    case 'ansi'
        
        R=zeros(1,length(F));             % response
        
        i=1;
        ref=log2(cross(1));
        while (F(i)<cross(1))
            R(i)=10^((ref-log2(F(i)))*(-12)/20); % -12
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
        
    otherwise
        error('Invalid type');
end

B=fir2(10000,F/fs*2, R);

if (bPlot)
    figure
    semilogx(F,20*log10(R),'-xr');
    title('Female filter');
    grid on;
    hold on
    [resp,Fresp]=freqz(B,1,[],fs);
    plot(Fresp,20*log10(abs(resp)), '-ob');
    
    %% Reference: SII
    [E f] = SpeechSptr('normal');
    plot(f,E-max(E),'g');
    
    legend('Target', 'Filter', 'ANSI S3.5 normal');
end

disp('')
