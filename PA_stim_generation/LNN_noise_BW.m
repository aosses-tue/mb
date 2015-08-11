function [LNN y] = LNN_noise_BW(N_iterations,varargin) % Fc,BW,SPL,dur,fs,Fmod,Mdept,dBFS)
% function [LNN y] = LNN_random_noise_BW(N_iterations,Fc,BW,SPL,dur,fs,Fmod,Mdept,dBFS)
% function [LNN y] = LNN_random_noise_BW(N_iterations,y,fs,Fc,BW)
% 
% 1. Description:
%       Creates one frame of AM-'running' noise at Fs, N.
%       The input parameters are centre frequency Fc and bandwidth BW
% 
%       For Gaussian-noise generation, see also: AM_random_noise_BW.m, AM_random_noise.m
% 
% 2.1 Example:
%       Finf = 950;
%       Fsup = 1050;
%       BW = Fsup - Finf;
%       fc = (Finf + Fsup)/2;
%       Mdept = 0; % if 0 then Fmod is not important
%       Fmod = 4;
%       dur = 0.5; % [s]
%       SPL = 70;
%       fs = 44100;
%       N_iterations = 10;
%       LNN = LNN_noise_BW(N_iterations,fc,BW,SPL,dur,fs,Fmod,Mdept);
%       sound(LNN,fs); 
% 
%       % If you want to store the output (Wav file). A reference plot will
%       % also be plotted as in Kohlrausch1997, Fig. 2
%       LNN_noise_BW(N_iterations,fc,BW,SPL,dur,fs,Fmod,Mdept); 
% 
% Programmed/adapted by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/08/2015
% Last update on: 06/08/2015
% Last use on   : 06/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    dBFS = 100;  
end
    
if nargin < 1
    N_iterations = 10;
end

if nargin > 2
    
    paramsName = {'Fc','BW','SPL','dur','fs','Fmod','Mdept','dBFS'}; % first parameter is common
    for i = 1:length(varargin)
        exp1 = sprintf('%s = varargin{i};',paramsName{i});
        eval(exp1);
    end
    
    if nargin < 8
        Mdept = 0; % no modulation as default
    end

    if nargin < 7
        Fmod  = 4;
    end

    if nargin < 6
        fs  = 40960;
    end

    if nargin < 5
        dur = 200e-3; % 200 ms, leading to N = 8192 at Fs = 40960 Hz
    end

    if nargin < 4
        SPL = 70;
    end

    if nargin < 3
        BW  = 100;
    end

    if nargin < 2
        Fc = 1000;
    end
    
    y = AM_random_noise_BW(Fc,BW,SPL,dur,fs,Fmod,Mdept,dBFS);
    
elseif nargin == 2
    y = varargin{1};
    fs = varargin{2};
    Fc = varargin{3};
    BW = varargin{4};
end

yenv = abs(hilbert(y));

fprintf('W=%.3f, V=%.1f (for no iteration)\n',mean(y.^4)/( mean(y.^2).^2 ),20*log10(std(yenv)/mean(yenv)));

t = ( 1:length(y) )/fs;
LNN = y;
Level_y = rmsdb(y);

finf = Fc - BW/2;
fsup = Fc + BW/2;

for i = 1:N_iterations
    
    env=abs(hilbert(LNN));
    LNN = LNN./env;
    
    LNNtmp = Set_Fourier_coeff_to_zero(LNN,fs,finf,fsup);
    
    LNN = setdbspl(LNNtmp,Level_y+dBFS); % set level to the same as input 'y'
    
    fprintf('W=%.3f, V=%.1f\n',mean(LNN.^4)/( mean(LNN.^2).^2 ), 20*log10(std(env)/mean(env)) );
    
end

LNNenv = env;
    
figure; 
subplot(3,1,1)
plot(t, y);
ylabel('Amplitude'); grid on

subplot(3,1,2)
plot(t,LNN,'r');
ylabel('Amplitude'); grid on

subplot(3,1,3)
plot(t, 20*log10(yenv), 'b',t,20*log10(LNNenv),'r');
ylabel('Envelope [dB]'); grid on

if nargout == 0
    
    sound(LNN,fs);
    filename = [Get_TUe_paths('outputs') sprintf('LNNnoise-fc-%.0f_BW-%.0f_Fmod-%.0f_Mdept-%.0f_SPL-%.0f',Fc,BW,Fmod,Mdept,SPL)];
    Wavwrite(LNN,fs,filename);
    
    filename = [Get_TUe_paths('outputs') sprintf('randomnoise-fc-%.0f_BW-%.0f_Fmod-%.0f_Mdept-%.0f_SPL-%.0f',Fc,BW,Fmod,Mdept,SPL)];
    Wavwrite(y,fs,filename);
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end