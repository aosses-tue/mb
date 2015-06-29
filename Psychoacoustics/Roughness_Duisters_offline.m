function [R dataOut out] = Roughness_Duisters_offline(insig, fs, N, options, CParams)
% function [R dataOut out] = Roughness_Duisters_offline(insig, fs, N, options, CParams)
%
% 1. Description:
%       Frame-based, off-line implementation of the roughness algorithm.
% 
%       Outputs:
%           out - has the same format as in the PsySound toolbox
% 
% author : Matt Flax <flatmax @ http://www.flatmax.org> : Matt Flax is flatmax
% March 2006 : For the psySoundPro project
%
% revised : Farhan Rizwi, July 2007
%           Reformatted, copied and vectorised code from InitAll
%           and Hweights into the function space below.  This
%           allows us to use nested functions effeciently.
%
% contact for the original source code :
%           http://home.tm.tue.nl/dhermes/
%
% 2. Stand-alone example:
%       2.1 Unix-based example:
%           [insig fs] = Wavread('~/Documenten/MATLAB/outputs/tmp-cal/ref_rough.wav'); % Unix-based
%           [out outPsy] = Roughness_Duisters_offline(insig,fs,8192);
% 
%       2.2 Multi-platform example, requires the ref_rough.wav file in the 
%           appropriate folder:
%           [insig fs] = Wavread([Get_TUe_paths('outputs') 'ref_rough.wav']); 
%           [out outPsy] = Roughness_Duisters_offline(insig,fs,8192);
% 
%       2.3 60-dB sine tone:
%           f = 1000;
%           fs = 44100;
%           lvl = 60;
%           insig = Create_sin(f,8192/fs,fs);
%           insig = setdbspl(insig,lvl);
%           [out outPsy] = Roughness_Duisters_offline(insig,fs,8192);
% 
%       2.4 Impulse response:
%           N = 8192;
%           insig  = [zeros(N/2-1,1); 1; zeros(N/2,1)];
%           fs = 44100;
%           [out outPsy] = Roughness_Duisters_offline(insig,fs,N);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 27/05/2015
% Last update on: 29/06/2015 % Update this date manually
% Last use on   : 29/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% switches(1) - Outer-, Middle- ear transmission
%       -  1: Van Immerseel
%       -* 2: Pfueger. This is the default
%
% switches(2) - Defines the auditory filterbank
%       -* 1: Gammatone, this is the default
%       -  2: Gamma chirp, not enabled yet
% 
% switches(3) - Defines an adaptation model
%       -* 2: Dau et al. 
%       - 10: same as 2, but with limitation of 10
%       - 12: same as 2, but gzi is determined as in Aures
% 
% switches(4) - Defines the resolution of the auditory filterbank + Roughness model to be used
%       -  1
%       -* 3: Aures

if nargin < 4
    options = [];
end

options = ef(options,'switches',[2 1 2 3]);

switches = options.switches;

% switches(1) = 2;
% switches(2) = 1; % or 2
% switches(3) = 2; % 2 or 10
% switches(4) = 3; % or 1, 3

if nargin < 3
    N = 8192;
end

if nargin < 5
    CParams = [];
    CParams.HopSize = N/2;
else
    CParams = Ensure_field(CParams,'HopSize',N/2); % 4096
end

options     = ef(options,'nSkipStart',0);
nSkipStart  = options.nSkipStart; % to correct time series, in case an analysis frame has been excluded
 
N_hop   = CParams.HopSize;
insig   = insig( nSkipStart*N_hop+1:end );
 

try
    fir_hweight; % just to check whether all folders are in the MATLAB path
catch
    warning('Temporal adjustment to add proper folders')
    misc.path = [Get_TUe_paths('MATLAB') 'Psychoacoustics' delim 'Roughness_Duisters' delim];
    Add_paths(misc);
end

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%

% BEGIN Hweights:
[numH7, denH7, numH14, denH14, numH30, denH30, numH36, denH36, numH66, denH66] = fir_hweight(fs);

% END Hweights 
% %%%%%%%%%%%%%%%%
 
insig = transpose(insig); % Transposing, just to follow Ronnies' nomenclature

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outer and middle-ear filtering
if     switches(1) == 1
% outer and middle ear filter as used by Van Immerseel and Martens
    disp('filtering outer and middle ear...');
    wr=2*pi* 4e3;
    num = wr^2;
    den = [1 0.33*wr wr^2];
    [num2, den2] = bilinear(num, den, fs);
    insig = filter(num2, den2, insig);
    
elseif switches(1) == 2
    disp('filtering outer and middle ear...');
    % Outer and middle ear combined bandpass filter
    % (Pflueger, Hoeldrich, Riedler, Sep 1997)
    % Low-pass component:
    b = 0.109*[1 1];
    a = [1 -2.5359 3.9295 -4.7532 4.7251 -3.5548 2.139 -0.9879 0.2836];
    % High-pass component
    d = [1 -2 1];
    c = [1 -2*0.95 0.95^2];
    num2 = conv(b, d);
    den2 = conv(a, c);
    insig = filter(num2, den2, insig);
    if fs ~= 48000
        warning('Pflueger''s outer-ear discrete approximation validated only for 48 kHz')
    end
else
    disp('no outer- and middle-ear filtering');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading defaults:

% auditory filterbanks
if     switches(4) == 1
    defin4.import={'auditoryfilterbank_Duisters1'}; % ERB 1:  1:  40
elseif switches(4) == 3 % This one was the default apparently
    defin4.import={'auditoryfilterbank_Duisters3'}; % ERB 0.5:0.5:38
end
[flags4,keyvals4]  = ltfatarghelper({'flow','fhigh'},defin4,{});
 
if     switches(3) == 2 | switches(3) == 12
    defin3.import           = {'ihcenvelope','adaptloop'}; 
    defin3.importdefaults   = {'ihc_breebaart','adt_dau1996'};  % ihc_breebaart (has a cuttoff of 770 Hz)
                                                                % adt_dau1996 it does not have limitation
elseif switches(3) == 10
    defin3.import           = {'ihcenvelope','adaptloop'}; 
    defin3.importdefaults   = {'ihc_breebaart','adt_dau'};  % ihc_breebaart (has a cuttoff of 770 Hz)
                                                            % adt_dau, lim = 10
end

[flags3,keyvals3]  = ltfatarghelper({'minlvl'},defin3,{});

%% Stage 1, BEGIN: RoughBody

% auditory filterbanks
if switches(2) == 1
    [inoutsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags4,keyvals4);
    erbr = freqtoaud(fc,'erb');
    Nch  = length(erbr);
end

% inner hair cells + adaptation models
if switches(3) == 2 | switches(3) == 10 | switches(3) == 12
    inoutsig        = ihcenvelope(inoutsig,fs,'argimport',flags3,keyvals3);
    inoutsig        =   adaptloop(inoutsig,fs,'argimport',flags3,keyvals3);
end
inoutsig = squeeze( inoutsig );

m_blocks = floor(length(insig)/N_hop) - 1;
ri       = zeros(m_blocks,Nch); % Memory allocation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_j = 1:m_blocks
 
    Ni = (idx_j - 1) * N_hop + 1;
    Nf = Ni+N-1;
    
    tn(idx_j) = Ni-1; % sample number to determine time
    
    if switches(4) == 3
        
        if idx_j == 1
            disp('extracting roughness...');
            % adapted Aures model
        end
        exc = transpose( inoutsig );
        
        % calculation of the DC values
        etmp = abs(exc);
        s0 = kron(ones(1, length(exc)), mean(etmp, 2));
        excd = etmp - s0;
        
        % Weighting filtering
        sBP(1:11,:)  = filter( numH7,  denH7, excd( 1:11,Ni:Nf), [], 2);
        sBP(12:28,:) = filter(numH14, denH14, excd(12:28,Ni:Nf), [], 2);
        sBP(29:36,:) = filter(numH30, denH30, excd(29:36,Ni:Nf), [], 2);
        sBP(37:65,:) = filter(numH36, denH36, excd(37:65,Ni:Nf), [], 2);
        sBP(66:76,:) = filter(numH66, denH66, excd(66:76,Ni:Nf), [], 2);
        
        sBPrms = rmsDik(sBP);
        rexc = rmsDik(exc);
        maxi = max(rexc);
        if maxi > 0
            calib = rexc / maxi;
        else
            calib = 0;
        end
        % modulation depth estimation
        for k = 1:Nch,
            % calibration factor
            if s0(k) > 0
                mdepth(k) = sBPrms(k) / s0(k);
                mdepth(k) = mdepth(k) * calib(k);
            else
                mdepth(k) = 0;
            end
            % calculation of the shifted cross correlation factor
            if k < Nch - 1,
                amount = 0.003 * fs;
                ki(k) = shiftcov(sBP(k,:), sBP(k+2,:), amount);
            end
        end
    end
    
    % definition of gzi
    if switches(2) == 1
        if switches(3) == 1 % gammatone,
            fc4gz   = [0 125  250 500  1000 2000 4000 8000 16000];
            Rmax    = [0 0.35 0.8 0.99 1    0.75 0.57 0.53 0.42]; % Van Immerseel & Martens
        elseif switches(3) == 2 
            fc4gz   = [0 125  250 500  1000 2000 4000 8000 16000];
            Rmax    = [0 0.3 1 1 1 0.64 0.49 0.51 0.45];
        elseif switches(3) == 3
            fc4gz   = [0 125  250 500  1000 2000 4000 8000 16000];
            Rmax    = [0 0.35 0.8 0.9 1 0.65 0.47 0.43 0.32]; % gammatone, Meddis
        elseif switches(3) == 12 | switches(3) == 10 
            gzitmp  = Get_psyparams('gr-ERB');
            fc4gz   = audtofreq( gzitmp(1,:),'erb' );
            Rmax    = gzitmp(2,:); 
        end
    end
    
    ERBrate = freqtoaud(fc4gz,'erb'); %2 * 21.4 .* log10(4.37 * fc / 1000 + 1);
    gzi     = interp1(ERBrate, Rmax, erbr, 'cubic');
    
    % figure;
    % plot(eb  ,gzzi(2,:)), hold on
    % plot(erbr,gzi,'r'); grid on
    % xlabel('Frequency [ERB]')
    
    % calculate specific roughness ri
    ri(idx_j,1:7)     = (gzi(1:7)     .* mdepth(1:7)     .* ki(1:7)).^2;
    ri(idx_j,8:Nch-3) = (gzi(8:Nch-3) .* mdepth(8:Nch-3) .* ki(6:Nch-5) .* ki(8:Nch-3)).^2;
    ri(idx_j,Nch-2)   = (gzi(Nch-2)    * mdepth(Nch-2)    * ki(Nch-4))^2;
    ri(idx_j,Nch-1)   = (gzi(Nch-1)    * mdepth(Nch-1)    * ki(Nch-3))^2;
    ri(idx_j,Nch)     = (gzi(Nch)      * mdepth(Nch)      * ki(Nch-2))^2;
    R(idx_j)          = sum(ri(idx_j,1:Nch));
     
end
 
% Create a cell array to return
dataOut{1}  = R;
dataOut{2}  = ri;
% out.t       = transpose(tn/fs);
 
nParam      = 1;

if nargout >= 3
    out.Data1   = transpose(R);
    out.name{nParam} = 'Roughness';
    out.param{nParam} = strrep( lower( out.name{nParam} ),' ','-');

    nParam      = 2;
    out.Data2   = ri;
    out.name{nParam} = 'Specific roughness';
    out.param{nParam} = strrep( lower( out.name{nParam} ),' ','-');
    out.t_ri    = transpose(tn/fs); 
    out.fi      = fc;
    out.erbi    = erbr;
    out.exc     = exc;
    out.nAnalyser = 15;
end

end
