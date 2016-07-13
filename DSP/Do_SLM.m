function [outsig_dB, outsig] = Do_SLM(insig,fs,weight_freq,weight_time)
% function [outsig_dB, outsig] = Do_SLM(insig,fs,weight_freq,weight_time)
%
% 1. Description:
%
% 2. Stand-alone example:
%   fs = 44100;
%   dur = 1; % s
%   insig = Create_sin(1000,dur,fs);
%   Do_SLM(insig,fs,'A','f');
% 
%   [insig, fs] = Wavread('D:\Databases\dir03-Speech\dutch\LISTman\jwz551.wav');
%   Do_SLM(insig,fs,'Z','f');
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/07/2016
% Last update on: 12/07/2016 
% Last use on   : 12/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    weight_freq = 'A';
end

if nargin < 4
    weight_time = 'f';
end

N = 8192*2;
fs = 44100;
weightings = il_gen_weighting_Filters(fs,N,weight_freq);
a = weightings.a;
b = weightings.b;

dBFS = 100;
curr_lvl = rmsdb(insig)+dBFS;
insig = setdbspl(insig,curr_lvl,'dboffset',94);

gain_dB = dBFS-94; % already applied

outsig = filter(b,a,insig);
outsig = il_integrator(outsig,fs,weight_time);

outsig_dB = 20*log10(abs(outsig/2e-5));

outsig = From_dB(-gain_dB)*outsig;

if nargout == 0
    figure;
    plot( outsig_dB ); grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function weightings = il_gen_weighting_Filters(fs,N,weightingType)
% function weightings = il_gen_weighting_Filters(fs,N,weightingType)
% 
% Generates the weighting filters

% the frequency of each Fourier bin
f = (0:fs/2/(N-1):fs/2);

% take a perceptual sampling of the frequency domain - for design purposes 
eMin   = freqtoaud(min(f)); % ERB min
eMax   = freqtoaud(max(f)); % ERB max
eScale = eMin:(eMax-eMin)/(length(f)-1):eMax; % ERB scale
fScale = audtofreq(eScale); % frequencies sample according to a linear ERB scale

fLinear = f; % save the linear frequency scale
f = fScale;  % switch the reference frequencies to be f

s = i*2*pi*f; % set up the s-plane variable

% determine the weighting filter frequency responses convienient to 
% accurately set the desired filter orders (n,m)
  
switch weightingType
    case 'A' % A-weighting filter
        K = 7.39705e9;
        freqResp = K*s.^4./((s+129.4).^2 .* (s+676.7).*(s+4636).*(s+76655).^2);
    
        n = 4; % at most we need a 4th order filter
        m = n; 
  
        zrs =  [0; 0; 0; 0];
        pls = -[129.4; 129.4; 676.7; 4636; 76655; 76655];
   
    case 'B' % B-weighting filter
        K = 5.99185e9;
        freqResp = K*s.^3./((s+129.4).^2 .* (s+995.9).*(s+76655).^2);

        n = 3; % at most we need a 4'th order filter
        m = 4; 
   
        zrs =  [0; 0; 0];
        pls = -[129.4; 129.4; 995.9; 76655; 76655];

    case 'C' % C-weighting filter
        K = 5.91797e9;
        freqResp = K*s.^2./((s+129.4).^2 .*(s+76655).^2);

        n = 2; % at most we need a 4'th order filter
        m = 4; 

        zrs =  [0; 0];
        pls = -[129.4; 129.4; 76655; 76655];

    case 'D' % D-weighting filter
        K = 91104.32;
        freqResp = K*s.*(s.^2+6532*s+4.0975e7)./((s+1776.3) .*(s+7288.5).*(s.^2+21514*s+3.8836e8));

        n = 3; % at most we need a 4'th order filter
        m = 4; 
   
        zrs = [0; roots([1 6532 4.0975e7])];
        pls = [-1776.3; -7288.5; roots([1 21514 3.8836e8])];
    
	 case 'R'
		% Filter weightings from [1], pg 12
		% These are defined for 48k
		b = [1 -2 1];
		a = [1 -1.99004745483398 0.99007225036621];

		% Use direct substituition of the definition of the z-transform
		% (z=exp(s*T)) to recalculate coeffecients for a different sampling
		% rate
		% Note: This could be another option for pre-filtering

		if Fs ~= 48e3;
            poles = roots(a);
  
            % Make polynomial after fixing up the roots
            % 
            % z = exp(s*T) --> s = ln(z)/T
            %
            % s = ln(z1)/T1 = ln(z2)/T2  -->  z2 = exp(ln(z1)*T2/T1)
            %
            a = poly(exp(log(poles)*48e3/Fs));

            % Note that the two zeros at 1 remain there.
            % Note also, that the negligible high frequency gain adjustment
            % is ignored.
		end

		weightings.a = a;
        weightings.b = b;
    
        % ... and we're done
        return		
		
    case 'Z' % un-weighted
        weightings.a = 1;
        weightings.b = 1;

        % ... and we're done
        return
    
    otherwise % unknown request
        error(['weightingType=''' weightingType ''' is unknown. Options are ''A'', ''B'', ''C'' or ''D'''])
end
  
m = m+1;
n = n*2;
m = m*2;
  
% the total frequency response
totalResp = freqResp;

% generate the filter
if (2*f ~= fs) % correct small frequency error on the last fourier sample.
    f(end) = fs/2;
end
  
% Use the bilinear transformation to discretize the above
% transfer function.

warnState = warning('off', 'MATLAB:nearlySingularMatrix');

[Zd, Pd, Kd] = bilinear(zrs, pls, K, fs);
[b, a] = zp2tf(Zd, Pd, Kd);
warning(warnState);

% Assign
weightings.a = a;
weightings.b = b;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_integrator(insig,fs, weightingType)
% function outsig = il_integrator(insig,fs, weightingType)
% 
% INTEGRATOR 
% 
% Generates and applies the integration filter to the input data
%
% For the theory behind this code, please refer to the document 'psy-sound.pdf'
%
% Author : Matt R. Flax <flatmax>
%          Jan. 2007 for the psysound project.
%
% Revised : Farhan Rizwi
%           Mostly cleanup of redundant code
%           July 2007 for the psysound project.
%
% input :
%         dataIn     - data vector
%         fastOrSlow - RC time constant is 'f' fast (125 ms) or 's'
%                      slow (1 s)
%         Fs         - sample rate of the data
%     
% Time constant for the leaky integrator is tau. This is basically
% a low-pass filter with a bandwidth of 1/tau.
%
% This yields the following transfer function :
%                      1
%        H(s) =  ---------------
%                tau s  +   1
%
% This can also be used with an rms integrator - 
% make fastorslow = 'r0.2' for an rms integrator with a window size of 0.2
% seconds. For a leaky integrator or 0.2 seconds 'l0.2' will work. 

% Filter coeffecients

switch weightingType
    case 'f' % fast leak - time constant = 125 ms
        tau = 125e-3;
    case 's' % slow leak - time constant = 1 s
        tau = 1;
    case 'i'
        tau = 35e-3; % impulse
    case 'p'
        tau = 50e-6;	
    case 'l'
        tau = str2num(fastOrSlow(2:end));	
    case 'r'
        tau = str2num(fastOrSlow(2:end));	
    otherwise
        error(['integrator: unknown leak case ' char(fastOrSlow)]);
end

% State vector
Z = [];

% Exponential term
E = exp(-1/(tau*fs));

% Filter numerator - with gain adjustment
b = 1 - E;

% Filter denominator
a = [1 -E];

% State vector
Z = [];

% Create run function handle
% Use filter to perform the integration
[outsig, Z] = filter(b, a, abs(insig), Z, 1);
