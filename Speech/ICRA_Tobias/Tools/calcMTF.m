function [mtf,cfModHz,cfHz] = calcMTF(input,fsHz,cfModHz,method,nSubbands)
%calcMTF   Compute the modulation transfer function (MTF). 
%   The envelope of the input signal is analyzed by a bank of modulation
%   filters [1,2]. The corresponding modulation analysis is either
%   implemented in the FFT domain [1,2] or using time-domain filters with a
%   constant Q-value [3,4]. The MTF can be computed based on the broadband
%   signal, or for a predefined number of subbands.       
% 
%USAGE
%   [mtf,cfModHz,cfHz] = calcMTF(input,fsHz)
%   [mtf,cfModHz,cfHz] = calcMTF(input,fsHz,cfModHz,method,nSubbands)
%
%INPUT ARGUMENTS
%       input : input signal [nSamples x 1]
%        fsHz : sampling frequency in Hz
%     cfModHz : modulation filter center frequencies in Hertz
%               (default, fRangeMod = [0.5 32])
%      method : string specifying method and resolution of modulation 
%               analysis '<method>_<resolution>' (default, method = 'fft_1') 
%               <method> = fft : compute FFT-based modulation spectrum and
%                                integrate the energy of FFT bins with
%                                equal weight across the frequency range
%                                corresponding to the individual modulation
%                                filters.   
% 
%                                'fft_<resolution>' controls the overlap
%                                between adjacent modulation filters  
%                                (e.g. <resolution> = 1 for 1-octave or 3
%                                for 1/3-octave modulation analysis)     
% 
%                          filter : measure the RMS at the output of a bank
%                                   of time-domain modulation filters.  
% 
%                                   'filter_<resolution>' controls the
%                                   quality factor of the modulation
%                                   filters (e.g. <resolution> = 1 for a
%                                   constant Q-value of 1 according to
%                                   [3,4]). A similar resolution compared
%                                   to 'fft_1' is obtained with filter when
%                                   the setting 'filter_1.4142' is used.       
% 
%   nSubbands : number of MEL-scaled subbands (default, nSubbands = 1)
% 
%OUTPUT ARGUMENTS
%         mtf : MTF for each subband and modulation filter [nSubbands x nFiltersMod]
%     cfModHz : modulation filter center frequencies [nFiltersMod x 1]
%        cfHz : subband filter center frequencies [nSubbands x 1]
% 
%   calcMTF(...) plots the MTF in a new figure.
% 
%   See also calcLTAS, calcLTAS_OCT, calcSTFT and calcPowerSpec.
% 
%REFERENCES
%   [1] W. A. Dreschler, H. Verschuure, C. Ludvigsen and S. Westermann,
%       "ICRA Noises: Artificial Noise Signals with Speech-like Spectral
%       and Temporal Properties for Hearing Instrument Assessment",
%       International Journal of Audiology, 40(3), pp.148-157, 2001. 
% 
%   [2] T. Houtgast and H. J. M. Steeneken, "The modulation transfer
%       function in room acoustics as a predictor of speech
%       intelligibility", Acustica, 28, pp.66-73, 1973.  
% 
%   [3] S. D. Ewert and T. Dau, "Characterizing frequency selectivity for
%       envelope fluctuations", The Journal of the Acoustical Society of
%       America, 108(3), pp.1181-1196, 2000.  
% 
%   [4] T. May and T. Dau, "Computational speech segregation based on an
%       auditory-inspired modulation analysis", The Journal of the
%       Acoustical Society of America, 136(6), pp.3350-3359, 2014.   

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/09
%   v.0.2   2015/02/20 added subband processing
%   v.0.3   2015/04/02 added FFT-based MTF calculation
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(cfModHz);   cfModHz   = [0.5 1 2 4 8 16 32]; end
if nargin < 4 || isempty(method);    method    = 'fft_1';             end % fft_3 for one-third octave band
if nargin < 5 || isempty(nSubbands); nSubbands = 1;                   end

% Check if input is mono
if min(size(input)) > 1
    error('Input must be mono.')
end

% Detect underscore 
idxSep = strfind(method,'_');

% Determine modulation analysis method
if strfind(method,'fft')
    if isempty(idxSep)
        nBandsPerOct = 1;
    else
        nBandsPerOct = str2double(method(idxSep(1)+1:end));
        if ~isfinite(nBandsPerOct);
            error('Wrong usage of <method>_<resolution>.')
        end
    end
    bUseFFT = true;
elseif strfind(method,'filter')
    if isempty(idxSep)
        qFactor = sqrt(2);
    else
        qFactor = str2double(method(idxSep(1)+1:end));
        if ~isfinite(qFactor);
            error('Wrong usage of <method>_<resolution>.')
        end
    end
    bUseFFT = false;
    bUseLP  = false; % Use band-pass filters only
    bUseHP  = false; % Use band-pass filters only
else
    error('MTF analysis method is not supported!')
end


%% NORMALIZE INPUT
% 
% Scale input to have zero mean and unit variance
input = normalizeData(input,'meanvar');

%% DECOMPOSE INPUT INTO SUBBAND SIGNALS
% 
% Compute subbands
if nSubbands > 1
    % Create 4th order MEL-scaled filterbank (first filter is a low-pass)
    MFB = createFB_MEL(fsHz,[50 fsHz/2],nSubbands,true,4);
        
    % Subband center frequencies
    cfHz = MFB.cfHz;
    
    % Apply filterbank
    subbands = applyFB(MFB,input);
else
    % Use broadband signal
    subbands = input; cfHz = [];
end


%% PERFORM ENVELOPE EXTRACTION
% 
% Perform intensity envelope extraction according to [2]
subbands = abs(subbands).^2;

% % Low-pass filter with cut-off 100 Hz and attenuation of 48 dB/oct
% [b,a] = Create_Butter_att_oct(100,fsHz,48,'low');
% subbands = filter(b,a,subbands);

%% DOWN-SAMPLE INPUT
% 
% Reference sampling frequency of the envelope signal after decimation
% (default => 1200 Hz). If higher modulation frequencies are requested,
% this reference will be changed accordingly.
fsHzEnvelope = max(1200,4 * round(max(cfModHz)));

% Resample absolute values of input
inputD = abs(resample(subbands,fsHzEnvelope,fsHz));


%% PERFORM MODULATION ANALYSIS
% 
% Number of modulation filters
nFiltersMod = numel(cfModHz);

% Allocate memory
pspecMod = zeros(nSubbands,nFiltersMod);
pspecDC  = zeros(nSubbands,1);


if bUseFFT
    % Determine FFT resolution
    nfft = pow2(nextpow2(size(inputD,1)));
    
    % Frequency vector from 0 to fsHz
    freqHz = (0:(nfft/2))'/(nfft/2)/2 * fsHzEnvelope;
        
     % Find lower and higher 3dB edge freqencies
    fModLowHz  = cfModHz .* 2^(-1 / nBandsPerOct / 2);
    fModHighHz = cfModHz .* 2^( 1 / nBandsPerOct / 2);

    % Bandwidth in Hertz
    bwHz = fModHighHz - fModLowHz;
    
    % Q-factor
    qFactor = cfModHz ./ bwHz; %#ok
    
    % Loop over the number of subbands
    for ii = 1 : nSubbands
        
        % Reduce discontinuities by windowing
        subband = fade(inputD(:,ii),fsHzEnvelope,10E-3,@hann); % 10-ms fade
        
        % FFT analysis
        spec = abs(time2freq(subband,nfft))/nfft;
        
        % Take positive frequencies times two 
        % (to reflect energy of negative frequencies)
        if rem(nfft,2)
            % Single-sided power spectrum is odd, only DC is unique
            spec(2:end,:) = spec(2:end,:) * 2;
        else
            % Single-sided power spectrum is even, do not double nyquist 
            spec(2:end-1,:) = spec(2:end-1,:) * 2;
        end
        
        % Loop over number of modulation filter
        for mm = 1 : nFiltersMod
            
            % Find frequency indices for mm-th modulation filter
            idxMod = fModLowHz(mm) < freqHz & freqHz < fModHighHz(mm);
            
            % Check if modulation filter is represented
            if sum(idxMod) == 0
                error(['No FFT bins detected for modulation filter '  ,...
                    'centered at ',num2str(cfModHz(mm)),' Hz. '       ,...
                    'Either increase the center frequency of lowest ' ,...
                    'modulation filter or increase the signal duration.'])
            end
        
            % Measure the root sum of squared values at the output of each
            % modulation filter
            pspecMod(ii,mm) = calcRSS(spec(idxMod),1);
        end
        
        % Compute DC component of the envelope signal
        pspecDC(ii) = mean(subband);
    end
else
    % Implement modulation filterbank
    ModFB = createFB_MOD(fsHzEnvelope,cfModHz,qFactor,bUseLP,bUseHP);
    
    % Loop over the number of subbands
    for ii = 1 : nSubbands
        
        % Reduce discontinuities by windowing
        subband = fade(inputD(:,ii),fsHzEnvelope,10E-3,@hann);
        
        % Apply modulation filterbank
        out = applyFB(ModFB,subband);
        
        % Measure the RMS at the output of each modulation filter
        pspecMod(ii,:) = calcRMS(out,1);
        
        % Compute DC component of the envelope signal
        pspecDC(ii) = mean(subband);
    end
end


%% COMPUTE MODULATION TRANSFER FUNCTION (MTF)
% 
% Modulation transfer function
mtf = bsxfun(@rdivide,pspecMod,pspecDC);


%% SHOW RESULT
% 
% Plot results
if nargout == 0
    if nSubbands == 1
        figure;
        plot(mtf,'.-');
        grid on;
        set(gca,'xtick',1:numel(cfModHz),'xticklabel',num2str(cfModHz','%1.1f'))
        title('Modulation transfer function')
        xlabel('Modulation frequency (Hz)')
        ylabel('Modulation transfer function (MTF)')
    else
        figure;
        contourf(mtf);
        axis xy
        colorbar;
        set(gca,'xtick',1:nFiltersMod,'xticklabel',num2str(cfModHz','%1.1f'))
        set(gca,'ytick',1:nSubbands,'yticklabel',num2str(cfHz','%1.0f'))
        title('2D Modulation transfer function')
        xlabel('Modulation frequency (Hz)')
        ylabel('Center frequency (Hz)')
    end
end


