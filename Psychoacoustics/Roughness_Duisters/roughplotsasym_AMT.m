function roughplotsasym_AMT(fs)
% function roughplotsasym_AMT(fs)
%
% 1. Description:
%       Same as roughplotsasym.m but using th AMT-compatible Roughness 
%       implementation.
% 
% This script generates roughness plots for ramped and damped signals
% as a function of modulation depth and centre frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    fs = 44100;
end

NDuisters       = 8192*5;
CParams.HopSize = NDuisters/2;
[insigref fss]  = Wavread([Get_TUe_paths('outputs') 'ref_rough.wav']); 
insigref        = resample(insigref,fss,fs);
R01             = Roughness_Duisters_offline(insigref(1:8192*5),fs,NDuisters,[],CParams);

bDebug = 0;
bDoRoughnessDH  = 0;
bDoRoughnessRD  = 1;

exp1 = sprintf('roughinit-at-%.0f-Hz',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switches(1) = 1; % outer and middle-ear filtering (1 = Van Immerseel and Martens; 2 = Pflueger et al)
% switches(2) = 1; % (1 = Gammatone; 2 = Gammachirp)
% switches(3) = 2; % Adaptation model (1 = Van Immerseel; 2 = Dau et al; 3 = Meddis)
% switches(4) = 3; % (3: ERBres = 1, ERBmax = 40; 0: ERBres = 0.5, ERBmax = 38)
% switches(5) = 0;
% switches(6) = 0;
% switches(7) = 3; % (3 = Aures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N   = 8192;
t   = 0:1/fs:1 - 1/fs; % 1-s tones
m   = [0.4 0.6 0.8]; % modulation index
fc  = [1000 5000 10000]; % carrier frequency [Hz]
fmod  = 70; % Hz
dBSPL = 60;

% 1. Stimuli generation:
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

for i = 1:length(fc)
    for j = 1:length(m),
        % damped
        signal = (1 + m(j) * (Esaw(i,:) / max(Esaw(i,:)))) .* cos(2 * pi * fc(i) .* t - pi/2);
        win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
        signal = [signal(1:length(win)) .* win ... % window signal
                  signal(length(win)+1:length(signal)-length(win)) ...
                  signal(length(signal)-length(win)+1:length(signal)) .* fliplr(win)];
        insig1 = From_dB(dBSPL-100)*signal;
        insigname = ['insig' num2str(i) num2str(j) '_damped'];   
        exp1 = sprintf('%s=insig1;',insigname);
        eval(exp1);
        
        % ramped
        signal = (1 + m(j) * (Esaw2(i,:) / max(Esaw2(i,:)))) .* cos(2 * pi * fc(i) .* t + pi/2);
        win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
        signal = [signal(1:length(win)) .* win ... % window signal
                  signal(length(win)+1:length(signal)-length(win)) ...
                  signal(length(signal)-length(win)+1:length(signal)) .* fliplr(win)];
        insig2 = From_dB(dBSPL-100)*signal;
        insigname = ['insig' num2str(i) num2str(j) '_ramped'];   
        exp1 = sprintf('%s=insig2;',insigname);
        eval(exp1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(fc),
    fprintf('\nfc = %4.0i\n', fc(i));
    for j = 1:length(m),
        % damped
        insigname = ['insig' num2str(i) num2str(j) '_damped'];
        exp2 = sprintf('insig1=%s;',insigname);
        eval(exp2);
        if bDoRoughnessDH
            Rtmp = Roughness_offline(insig1(1:N),fs,N,0);
            RDH(i,j) = Rtmp(1);
        end
        
        if bDoRoughnessRD
            [Rtmp1 dataOut out] = Roughness_Duisters_offline(insig1, fs, NDuisters,[],CParams);
            RRD(i,j) = Rtmp1(1);
            
            if bDebug
                figure(5)
                subplot(2,1,1)
                idx = 30;
                plot( t, out.exc(idx,:) ); hold on
                title('damped')
            
                figure(6)
                plot( out.erbi, dataOut{1,2} ); hold on
            end
        end
        
        % ramped
        insigname = ['insig' num2str(i) num2str(j) '_ramped'];
        exp2 = sprintf('insig2=%s;',insigname);
        eval(exp2);
        if bDoRoughnessDH
            Rtmp = Roughness_offline(insig2(1:N),fs,N,0);
            RDH2(i,j)= Rtmp(1);
        end
        if bDoRoughnessRD
            [Rtmp1 dataOut out] = Roughness_Duisters_offline(insig2, fs, NDuisters,[],CParams);
            RRD2(i,j) = Rtmp1(1);
            
            if bDebug
                figure(5)
                subplot(2,1,2)
                plot( t, out.exc(idx,:) )
                title('ramped')

                figure(6)
                plot( out.erbi, dataOut{1,2} ,'r')
                legend('damped','ramped')
            end
        end
    end
end

if bDoRoughnessRD
    % convert scale to asper
    RRD  = RRD  / R01;
    RRD2 = RRD2 / R01;
    % auto axis scaling
    minax = min(min(min(RRD)), min(min(RRD2)));
    maxax = max(max(max(RRD)), max(max(RRD2)));
    minax = minax - 0.1 * max(abs(minax), abs(maxax));
    maxax = maxax + 0.1 * max(abs(minax), abs(maxax));
    legtop = maxax - 0.1 * maxax;
    figure(11);
    clf;

    subplot(2,3,1), plot(m, RRD(1,:), 'k.-', m, RRD2(1,:), 'ko--');
    xlabel('modulation depth');
    ylabel('roughness (asper)');
    axis([0.3 0.9 minax maxax]);
    text(0.35, legtop, ['f_c = ' num2str(fc(1)/1000) ' kHz']);
    subplot(2,3,2), plot(m, RRD(2,:), 'k.-', m, RRD2(2,:), 'ko--');
    xlabel('modulation depth');
    axis([0.3 0.9 minax maxax]);
    text(0.35, legtop, ['f_c = ' num2str(fc(2)/1000) ' kHz']);
    subplot(2,3,3), plot(m, RRD(3,:), 'k.-', m, RRD2(3,:), 'ko--');
    xlabel('modulation depth');
    legend('R_d', 'R_r', 1);
    axis([0.3 0.9 minax maxax]);
    text(0.35, legtop, ['f_c = ' num2str(fc(3)/1000) ' kHz']);
    
    figure;
    subplot(2,1,1)
    plot(t,insig23_damped);
    ha = gca;

    subplot(2,1,2)
    plot(t,insig23_ramped);
    ha(end+1) = gca;

    linkaxes(ha,'xy');
    xlim([0 0.06])
end
