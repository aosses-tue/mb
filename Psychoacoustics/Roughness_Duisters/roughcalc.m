%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auditory roughness estimation using several modules
%
% Author: Ronnie Duisters               %
% Group: Human-Technology Interaction   %
% Department of Technology Management   %
% Eindhoven University of Technology    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% R = roughcalc(signal, L, fs, winlength, dN, alpha, switches)
%
% inputs:
% signal input signal
% L preferred input level in dB SPL (default = 60)
% fs sample frequency in Hz (default = 48000)
% winlength length of analysis window in s (default = 0.2)
% dN sample interval to be processed by the
% roughness extraction model
% (default = [round(0.3*fs) round(0.5*fs)])
% alpha roughness index (default = 1.6)
% switches array selecting the preferred models
% (default = [2213])
%
% outputs:
% R roughness value
%
% switches [a b c d]
% a = middle ear filter:
% 0 = no filter
% 1 = Van Immerseel and Martens
% 2 = Pflueger et al.
% b = filterbank:
% 1 = Gammatone filterbank
% 2 = Gammachirp filterbank
% c = hair cell model:

% 0 = no hair cell model
% 1 = Van Immerseel and Martens hair cell model
% 2 = Dau et al. adaptation model
% 3 = Meddis hair cell model
% d = Roughness model
% 1 = SIM model I
% 2 = SIM model II
% 3 = modified Aures model (Jourdes)

function R = roughcalc(signal, L, fs, winlength, dN, alpha, switches)

global numH7 denH7 numH14 denH14 numH30 denH30 numH36 denH36 numH66 denH66

if nargin < 7
    switches = [2213];
end

if nargin < 6
    alpha = 1.6;
end

if nargin < 5
    dN = [round(0.3*fs) round(0.5*fs)];
end

if nargin < 4
    winlength = 200e-3;
end

if nargin < 3
    fs = 48000;
end

if nargin < 2
    L = 60;
end

if nargin < 1
    help roughcalc;
    error('no input signal defined!');
end

winlength = winlength * fs;
wsize = round(min(length(signal), winlength));
signal = signal(:)'; % make signal a row-vector
signal = signal(1:wsize);
disp('adapting signal level...');
signal = AdaptLevel(signal, L, switches);

% outer and middle-ear filtering
if switches(1) == 1
% outer and middle ear filter as used by Van Immerseel and Martens
    disp('filtering outer and middle ear...');
    wr=2*pi* 4e3;
    num = wr^2;
    den = [1 0.33*wr wr^2];
    [num2, den2] = bilinear(num, den, fs);
    signal = filter(num2, den2, signal);
elseif switches(1) == 2
    disp('filtering outer and middle ear...');
    % Outer and middle ear combined bandpass filter
    % (Pflueger, Hoeldrich, Riedler, Sep 1997)
    % Highpass component
    b = 0.109*[1 1];
    a = [1 -2.5359 3.9295 -4.7532 4.7251 -3.5548 2.139 -0.9879 0.2836];
    % Lowpass component
    d = [1 -2 1];
    c = [1 -2*0.95 0.95^2];
    signal = filter(conv(b, d), conv(a, c), signal);
else
    disp('no outer- and middle-ear filtering');
end

% auditory filterbanks
if switches(4) == 3
    erbres = 0.5;
    erbmax = 38;
else
    erbres = 1;
    erbmax = 40;
end

erbr = erbres:erbres:erbmax;
fc = 229 * (10.^(erbr / 21.4) - 1);
Nch = length(erbr);
if switches(2) == 1
    % gammatone filterbank
    disp('performing gammatone filtering...');
    fcoefs = MakeERBFilters(fs, fc(length(fc):-1:1), 26); %MakeERBFilters(fs, fc(length(fc):-1:1), 26, 1.019);
    fcoefs = fcoefs(Nch:-1:1,:);
    exc = ERBFilterBank(signal, fcoefs); % replace it for AMT function
elseif switches(2) == 2
    % compressive gammachirp filterbank
    disp('performing gammachirp filtering...');
    n=4;
    b1 = 1.81;
    c1 = -2.96;
    b2 = 2.17;
    c2 = 2.20;
    fr1 = fc;
    fp1 = fr1 + c1 * b1 * (24.7 + 0.107939 * fr1) / n;
    frat = 0.573 + 0.0101 * L;
    fr2 = frat * fp1;
    fcoefs = MakeERBFilters(fs, fc(length(fc):-1:1), 26, b1);
    fcoefs = fcoefs(Nch:-1:1,:);
    exc = ERBFilterBank(signal, fcoefs);
    for i = 1:length(exc(:,1)),
        [bgc1, agc1] = MakeAsymCmpFilters(fs, fr1(i), n, b1, c1);
        exc(i,:) = filter(bgc1, agc1, exc(i,:));
        [bgc2, agc2] = MakeAsymCmpFilters(fs, fr2(i), n, b2, c2);
        exc(i,:) = filter(bgc2, agc2, exc(i,:));
    end
    if level == 80, % low-frequency channels become unstable for the 80 dB condition
        for v = 1:2,
            exc(v,:) = 0 * exc(v,:);
        end
    end
end

% adaptation models
if switches(3) == 1
    % hair cell model by Van Immerseel and Martens
    disp('applying vI&M hair cell model...');
    y0 = 0.4472;
    exc = max(exc + y0, 0);
    fspont = 0.05e3;
    fsat = 0.15e3;
    r = 0.86;
    tau1 = 8e-3;
    tau2 = 40e-3;
    num = [tau1-r*tau1+r*tau2 1];
    den = [tau1*tau2 tau1+tau2 1];
    [num2, den2] = bilinear(num, den, fs);
    for i = 1:length(exc(:,1)),
        q(i,:) = filter(num2, den2, exc(i,:));
        f(i,:) = (fsat * exc(i,:)) ./ (((sqrt(fsat/fspont) - 1) * ...
        sqrt(y0)) + sqrt(q(i,:))).^2;
    end
    [bn, an] = butter(3, 600 / (fs / 2));
    exc=f;
    exc = filter(bn, an, exc);
    exc = max(exc-50, 0);
elseif switches(3) == 2
    % adaptation model by Dau et al.
    disp('applying Dau et al. adaptation model...');
    
    % Half-wave rectification:
    exc = max(exc, 0);
    [bn, an] = butter(5, 770 / (fs / 2));
    exc = filter(bn, an, exc, [], 2);
    
    CalFactor = -100;
    exc = From_dB(CalFactor) * exc;
    % for i = 1:Nch,
    %     f(i,:) = fadapt(exc(i,:)', fs)';
    % end
    % exc=f;
    % exc = max(exc, 0);
    
    limit = 0;
    minlvl = 1e-5;
    tau = [0.005 0.050 0.129 0.253 0.500];
    for i = 1:Nch
        f(i,:) = comp_adaptloop(exc(i,:)', fs, limit,minlvl,tau)';
        exc(i,:) = f(i,:);
    end
    
elseif switches(3) == 3
    % hair cell model by Meddis
    disp('applying Meddis hair cell model...');
    exc = max(exc, 0);
    exc = MeddisHairCell(exc, fs);
    exc = max(exc-65, 0);
else
    disp('no hair cell model');
end

% roughness extraction models
if switches(7) == 1
    % Leman roughness extraction model I
    disp('estimating roughness...');
    w = kron(ones(Nch, 1), hamming(wsize)');
    wexc = w .* exc(:,1:wsize);
    N = nextpow2(wsize);
    wsize = 2^N;
    fbin = fs / wsize;
    F = FilterWeights(erbmax/erbres, round(300/fbin), fc, 1, round(300/fbin), 1);
    if switches(2) == 0
        if switches(3) == 1 % gammatone,
            Rmax = [0 1.1 1.4 1.5 1 0.71 0.56 0.41 0.32]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 1.1 1.4 1.6 1 0.65 0.54 0.38 0.26]; % gammatone, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.9 1 1.2 1 0.65 0.5 0.38 0.26]; % gammatone, Meddis
        end
    elseif switches(2) == 2
        if switches(3) == 1 % gammachirp,
            Rmax = [0 0.3 1.3 1.4 1 0.74 0.73 0.42 0.38]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 0.3 0.94 1.4 1 0.7 0.6 0.38 0.31]; % gammachirp, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.84 0.94 1 1 0.82 0.64 0.44 0.35]; % gammachirp, Meddis
        end
    end
    fc = [0 125 250 500 1e3 2e3 4e3 8e3 16e3];
    ERBrate = 21.4 .* log10(4.37 * fc / 1000 + 1);
    gzi = interp1(ERBrate, Rmax, erbr, 'cubic');
    F = diag(gzi, 0) * F;
    l = length(F(1,:));
    sexc = fft(wexc, wsize, 2);
    sexc = sexc(:,1:l);
    DC = kron(sexc(:,1) / 2, ones(1, l));
    Bd=F.*sexc ./ DC;
    Rf = sum(abs(Bd) .^ alpha, 1);
    R = sum(Rf) / Nch;
    R = abs(R);
    Rf = abs(Rf);
elseif switches(7) == 2
    % Leman roughness extraction model II
    disp('estimating roughness...');
    w = kron(ones(Nch, 1), hamming(wsize)');
    wexc = w .* exc(:,1:wsize);
    N = nextpow2(wsize);
    wsize = 2^N;
    fbin = fs / wsize;
    F = FilterWeights(erbmax/erbres, round(300/fbin), fc, 1, round(300/fbin), 1);
    if switches(2) == 0
        if switches(3) == 1 % gammatone,
            Rmax = [0 1.1 1.4 1.5 1 0.71 0.56 0.41 0.32]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 1.1 1.4 1.6 1 0.65 0.54 0.38 0.26]; % gammatone, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.9 1 1.2 1 0.65 0.5 0.38 0.26]; % gammatone, Meddis
        end
    elseif switches(2) == 2
        if switches(3) == 1 % gammachirp,
            Rmax = [0 0.3 1.3 1.4 1 0.74 0.73 0.42 0.38]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 0.3 0.94 1.4 1 0.7 0.6 0.38 0.31]; % gammachirp, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.84 0.94 1 1 0.82 0.64 0.44 0.35]; % gammachirp, Meddis
        end
    end
    fc = [0 125 250 500 1e3 2e3 4e3 8e3 16e3];
    ERBrate = 21.4 .* log10(4.37 * fc / 1000 + 1);
    gzi = interp1(ERBrate, Rmax, erbr, 'cubic');
    F = diag(gzi, 0) * F;
    l = length(F(1,:));
    sexc = fft(wexc, wsize, 2);
    sexc = sexc(:,1:l);
    DC = sum(kron(sexc(:,1) / 2, ones(1, l)));
    Bd=F.*sexc;
    Bd2 = sum(Bd, 1) ./ DC;
    Rf = abs(Bd2) .^ alpha;
    R = sum(Rf) / Nch;
elseif switches(4) == 3
    % adapted Aures model
    exc = exc(:, dN(1):dN(2));
    disp('extracting roughness...');
    % calculation of the DC values
    etmp = abs(exc);
    s0 = kron(ones(1, length(exc)), mean(etmp, 2));
    excd = etmp - s0;
    % Weighting filtering
    sBP(1:11,:)  = filter(numH7 , denH7 , excd(1:11,:) , [], 2);
    sBP(12:28,:) = filter(numH14, denH14, excd(12:28,:), [], 2);
    sBP(29:36,:) = filter(numH30, denH30, excd(29:36,:), [], 2);
    sBP(37:65,:) = filter(numH36, denH36, excd(37:65,:), [], 2);
    sBP(66:76,:) = filter(numH66, denH66, excd(66:76,:), [], 2);
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
            mdepth(i) = 0;
        end
        % calculation of the shifted cross correlation factor
        if k < Nch - 1,
            amount = 0.003 * fs;
            ki(k) = shiftcov(sBP(k,:), sBP(k+2,:), amount);
        end
    end
    % definition of gzi
    if switches(2) == 1
        if switches(3) == 1 % gammatone,
            Rmax = [0 0.35 0.8 0.99 1 0.75 0.57 0.53 0.42]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 0.3 1 1 1 0.64 0.49 0.51 0.45]; % gammatone, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.35 0.8 0.9 1 0.65 0.47 0.43 0.32]; % gammatone, Meddis
        end
    elseif switches(2) == 2
        if switches(3) == 1 % gammachirp,
            Rmax = [0 0.35 0.8 0.99 1 0.88 0.64 0.68 0.65]; % Van Immerseel & Martens
        elseif switches(3) == 2
            Rmax = [0 0.05 0.15 0.99 1 0.73 0.61 0.77 0.6]; % gammachirp, Dau et al.
        elseif switches(3) == 3
            Rmax = [0 0.35 0.85 0.9 1 0.82 0.55 0.56 0.5]; % gammachirp, Meddis
        end
    end
    fc = [0 125 250 500 1e3 2e3 4e3 8e3 16e3];
    ERBrate = 2 * 21.4 .* log10(4.37 * fc / 1000 + 1);
    gzi = interp1(ERBrate, Rmax, erbr, 'cubic');
    % calculate specific roughness ri
    ri(1:7) = (gzi(1:7) .* mdepth(1:7) .* ki(1:7)).^2;
    ri(8:Nch-3) = (gzi(8:Nch-3) .* mdepth(8:Nch-3) .* ki(6:Nch-5) .* ki(8:Nch-3)).^2;
    ri(Nch-2) = (gzi(Nch-2) * mdepth(Nch-2) * ki(Nch-4))^2;
    ri(Nch-1) = (gzi(Nch-1) * mdepth(Nch-1) * ki(Nch-3))^2;
    ri(Nch) = (gzi(Nch) * mdepth(Nch) * ki(Nch-2))^2;
    R = sum(ri);
end
fprintf('\n');

end
