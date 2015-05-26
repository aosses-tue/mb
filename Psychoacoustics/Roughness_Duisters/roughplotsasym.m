% This script generates roughness plots for ramped and damped signals
% as a function of modulation depth and centre frequency.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
global numH7 denH7 numH14 denH14 numH30 denH30 numH36 denH36 numH66 denH66

fs = 48000;

exp1 = sprintf('roughinit-at-%.0f-Hz',fs);
exp2 = sprintf('load %s;',exp1);
try
    eval(exp2);
catch
    switches(1) = 1; % outer and middle-ear filtering (1 = Van Immerseel and Martens; 2 = Pflueger et al)
    switches(2) = 1; % (1 = Gammatone; 2 = Gammachirp)
    switches(3) = 2; % Adaptation model (1 = Van Immerseel; 2 = Dau et al; 3 = Meddis)
    switches(4) = 3; % (3: ERBres = 1, ERBmax = 40; 0: ERBres = 0.5, ERBmax = 38)
    switches(5) = 0;
    switches(6) = 0;
    switches(7) = 3; % (3 = Aures)
    N = 8192;
    if switches(4) == 3
        twin = 0.6; % N/fs;
        dN = [floor(0.0*fs) floor(0.0*fs)+N-1];
    else
        twin = 0.6;
        dN = [floor(0.3*fs) floor(0.5*fs)];
    end
    
    roughinit(fs, switches, twin, dN);
    movefile('roughinit.mat',[exp1 '.mat']);
    eval(exp2);
end

N   = 8192;
t   = 0:1/fs:0.6;
m   = [0.4 0.6 0.8];
fc  = [1000 5000 10000]; % Hz
fmod  = 70; % Hz
dBSPL = 60;

if switches(7) == 1
    alpha = 1.6;
else
    alpha = [];
end

% define observation interval

if switches(4) == 3
    twin = 0.6; % N/fs;
	dN = [floor(0.3*fs) floor(0.3*fs)+N-1];
else
    twin = 0.6;
    dN = [floor(0.3*fs) floor(0.5*fs)];
end
Esaw = zeros(3, length(t));
Esaw2 = zeros(3, length(t));
for i = 1:2,
    Esaw(1,:) = (1/i) * cos(2 * pi*i*fmod .* t - pi/2) + Esaw(1,:);
end
for i = 1:4,
    Esaw(2,:) = (1/i) * cos(2 * pi*i*fmod .* t - pi/2) + Esaw(2,:);
end
for i = 1:7,
    Esaw(3,:) = (1/i) * cos(2 * pi*i*fmod .* t - pi/2) + Esaw(3,:);
end
for i = 1:2,
    Esaw2(1,:) = (1/i) * cos(2 * pi*i*fmod .* t + pi/2) + Esaw2(1,:);
end
for i = 1:4,
    Esaw2(2,:) = (1/i) * cos(2 * pi*i*fmod .* t + pi/2) + Esaw2(2,:);
end
for i = 1:7,
    Esaw2(3,:) = (1/i) * cos(2 * pi*i*fmod .* t + pi/2) + Esaw2(3,:);
end

for i = 1:length(fc),
    fprintf('\nfc = %4.0i\n', fc(i));
    for j = 1:length(m),
        % damped
        signal = (1 + m(j) * (Esaw(i,:) / max(Esaw(i,:)))) .* cos(2 * pi * fc(i) .* t - pi/2);
        win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
        signal = [signal(1:length(win)) .* win ... % window signal
                  signal(length(win)+1:length(signal)-length(win)) ...
                  signal(length(signal)-length(win)+1:length(signal)) .* fliplr(win)];
    
        signal1 = From_dB(dBSPL-100)*signal;
        Rtmp = Roughness_offline(signal1,fs,N,0);
        RD(i,j) = Rtmp(1);
        
        R1(i,j) = roughcalc(signal, dBSPL, fs, twin, dN, alpha, switches);
        
        % ramped
        signal = (1 + m(j) * (Esaw2(i,:) / max(Esaw2(i,:)))) .* cos(2 * pi * fc(i) .* t + pi/2);
        win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
        signal = [signal(1:length(win)) .* win ... % window signal
        signal(length(win)+1:length(signal)-length(win)) ...
        signal(length(signal)-length(win)+1:length(signal)) .* fliplr(win)];
        
        signal2 = From_dB(dBSPL-100)*signal;
        Rtmp = Roughness_offline(signal2,fs,N,0);
        RD2(i,j)= Rtmp(1);
        R2(i,j) = roughcalc(signal, dBSPL, fs, twin, dN, alpha, switches);
        
        % figure; plot(signal1), hold on; plot(signal2,'r');
    end
end
% convert scale to asper
R1 = R1 / R01;
R2 = R2 / R01;
% auto axis scaling
minax = min(min(min(R1)), min(min(R2)));
maxax = max(max(max(R1)), max(max(R2)));
minax = minax - 0.1 * max(abs(minax), abs(maxax));
maxax = maxax + 0.1 * max(abs(minax), abs(maxax));
legtop = maxax - 0.1 * maxax;
figure(11);
clf;

subplot(231), plot(m, R1(1,:), 'k.-', m, R2(1,:), 'ko--');
xlabel('modulation depth');
ylabel('roughness (asper)');
axis([0.3 0.9 minax maxax]);
text(0.35, legtop, 'f_c = 2.5 kHz');
subplot(232), plot(m, R1(2,:), 'k.-', m, R2(2,:), 'ko--');
xlabel('modulation depth');
axis([0.3 0.9 minax maxax]);
text(0.35, legtop, 'f_c = 5 kHz');
subplot(233), plot(m, R1(3,:), 'k.-', m, R2(3,:), 'ko--');
xlabel('modulation depth');
legend('R_d', 'R_r', 1);
axis([0.3 0.9 minax maxax]);
text(0.35, legtop, 'f_c = 10 kHz');

disp('')