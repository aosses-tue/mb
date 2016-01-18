function [dataOut out] = FluctuationStrength_Garcia_offline(insig, fs, N, optsDebug)
% function [dataOut out] = FluctuationStrength_Garcia_offline(insig, fs, N, optsDebug)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
% 
%       Comments on this implementation: cross-correlation seems not to be 
%       working properly.
%
%       Adapted from Model/Helper/FluctuationStrength.m (first adaptation from Roughness model)
%       Adapted from Model/Helper/FluctuationStrength_debug.m (main file as in the thesis)
% 
%       Comments:
%           fs  - should be an input parameter
%           N   - should be an input parameter
%           Too difficult to change N!
%           gzi - where was it taken from?
% 
% 2. Stand-alone example:
%       r20141126_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Rodrigo Garcia L./Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 29/06/2015
% Last update on: 16/11/2015 
% Last use on   : 16/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    optsDebug = [];
end

optsDebug   = ef(optsDebug,'all',0);
optsDebug   = ef(optsDebug,'ki',0);

bDebug      = optsDebug.all;
debug_ki    = optsDebug.ki;

if optsDebug.all == 1 | optsDebug.ki == 1
    disp('Debug mode active. Only the first frame will be analysed');
    hFig = [];
end

% General fluctuation strength parameters
Ndataset = 0;
params  = FluctuationStrength_Garcia_getparams(N,fs,Ndataset);

if nargin < 3
    N   = params.N;
end

a0      = params.a0;
Barkno  = params.Barkno;
Cal     = params.Cal;
Chno    = params.Chno;
% MinBf   = params.MinBf; % Delete from _getparams
MinExcdB= params.MinExcdB;
N01     = params.N01;
p_g     = params.p_g;
p_m     = params.p_m;
p_k     = params.p_k;
qb      = params.qb;
zb      = params.zb;
gzi     = il_create_gzi(Chno);
Hweight = params.Hweight;
    
% Multiframe separation
overlap = 0.50;
b       = buffer(insig,N,N * overlap,'nodelay');
nFrames = size(b,2);
    
% Blackman window
window  = blackman(N);
ampCal  = From_dB(100) * 2 / (N * mean(window));
window  = ampCal * window';
    
FS      = zeros(1,nFrames);

for iFrame = 1:nFrames
    
    h0      = zeros(   1,Chno);
    hBPi    = zeros(Chno,   N);
    mdept	= zeros(   1,Chno);
    
    x = b(:,iFrame);

    SPL(iFrame) = rmsdb(x)+100;
    tempIn  = window .* x';
    
    %% Peripheral stage:
    
    % Transmission factor a0:
    FreqIn  = a0 .* fft(tempIn);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Critical-band filterbank - Terhardt:
    [etmp_td, etmp_fd, etmpExc_fd] = Terhardt_filterbank(FreqIn,fs,N,qb,MinExcdB,Barkno,N01,zb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear etmp_fd etmpExc_fd; % it liberates some memory
    
    etmp_td = abs(etmp_td);
    h0		= mean( etmp_td ,2); % mean along dimension 2
    
    for k = 1:Chno
        
        % Fei         = fft( etmp_td(k,:) - h0(k) );
        % hBPi(k,:)   = 2 * real(ifft(Fei.* Hweight'));
        hBPi(k,:) = 2*filter(Hweight,etmp_td(k,:)-h0(k));
        
    end
    
    hBPrms  = rms( hBPi ,'dim',2);
    
    idx = find(h0 > 0);
    mdept(idx) = hBPrms(idx) ./ h0(idx);
    
    idx = find(mdept > 1);
    mdept(idx) = 1;
    
    %%%
    warning('Temporal')
    Md = mdept - 0.1;
    idx = find(Md < 0);
    Md(idx) = 0;
    %%%
    
    ki = zeros(1,Chno - 2);
    fi = zeros(1,Chno);

    % Find cross-correlation coefficients
    for k=1:1:Chno-2
        cfac    = cov(hBPi(k,:),hBPi(k + 2,:));
        den     = diag(cfac);
        den     = Round( sqrt(den * den'), 10); % rounded to 10 decimals

        if den(2,1) > 0
            ki(k) = cfac(2,1) / den(2,1);
        elseif den(2,1) < 0
            ki(k) = 0;
        else
            ki(k) = 0;
        end
    end

    % Calculate specific fluctuation strength fi and total FS
    fi(iFrame,1) = gzi(1)^p_g * Md(1)^p_m * abs(ki(1))^p_k;
    fi(iFrame,2) = gzi(2)^p_g * Md(2)^p_m * abs(ki(2))^p_k;

    for k = 3:1:45
        fi(iFrame,k) = gzi(k)^p_g * Md(k)^p_m * abs( ki(k - 2) * ki(k) )^p_k;
    end

    fi(iFrame,46) = gzi(46)^p_g * Md(46)^p_m * abs(ki(44))^p_k;
    fi(iFrame,47) = gzi(47)^p_g * Md(47)^p_m * abs(ki(45))^p_k;

    FS(iFrame) = Cal * sum(fi(iFrame,:));
end

dataOut{1} = FS;
dataOut{2} = fi;
dataOut{3} = SPL;

nParam      = 1;
out.Data1   = transpose(FS);
output.name{nParam} = 'Fluctuation strength';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

function gzi = il_create_gzi(Chno);

    Chstep = 0.5;
    
    g0 = [
        0       0
        125     1
        250     1
        500     1
        1000    1
        1500    1
        2000    1
        3000    1
        4000    1
        6000    1
        8000    1
        16000   0
    ];
    
    gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
    gzi(isnan(gzi)) = 0;
    
    gzi = ones(1,Chno);