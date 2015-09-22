function r20150910_roughness_as_function_of
% function r20150910_roughness_as_function_of
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 10/09/2015
% Last update on: 10/09/2015 
% Last use on   : 10/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPart1 = 0;
bPart2 = 1;

N = 8192;
if bPart1
[insig fs] = Wavread([Get_TUe_paths('outputs') 'ref_rough.wav']);
insigm20 = gaindb(insig, -20);
insigp20 = gaindb(insig,  20);

[R(2) outPsy out2] = Roughness_offline(insig(1:N),fs,N);

[R(1) outPsy out1] = Roughness_offline(insigm20(1:N),fs,N);
[R(3) outPsy out3] = Roughness_offline(insigp20(1:N),fs,N);

figure;
subplot(3,1,1)
plot(out1.ri(1,:),'r'), hold on
plot(out1.mdept,'b-')
plot(out1.ki,'k')
plot(out1.gzi,'g--')
grid on;
title('test level at 40 dB SPL')

subplot(3,1,2)
plot(out2.ri(1,:),'r'), hold on
plot(out2.mdept,'b-')
plot(out2.ki,'k')
plot(out2.gzi,'g--')
grid on;
title('test level at 60 dB SPL')

subplot(3,1,3)
plot(out3.ri(1,:),'r'), hold on
plot(out3.mdept,'b-')
plot(out3.ki,'k')
plot(out3.gzi,'g--')
grid on;
title('test level at 80 dB SPL')

legend('r_i','m','cc','gzi')
% % Review: different values than Roughness_offline
% [Rd(2)] = Roughness_offline_debug(insig(1:N),fs,N);
% [Rd(1)] = Roughness_offline_debug(insigm20(1:N),fs,N);
% [Rd(3)] = Roughness_offline_debug(insigp20(1:N),fs,N);
end

SPL = 60;
if bPart2
fs = 44100;
fc = 1000;
dur = N/fs;
fmod = 70;
m = [0 0.5 1]*100;
outam = [];

figure;
for i = 1:3
    insigAM(:,i) = AM_sine(fc,dur,fs,fmod,m(i),SPL);
    [Ram(i) outPsy out] = Roughness_offline(insigAM(1:N,i),fs,N);
    
    subplot(3,1,i)
    plot(out.ri(1,:),'r'), hold on
    plot(out.mdept,'b-')
    plot(out.ki,'k')
    plot(out.gzi,'g--')
    grid on;
    title(sprintf('Modulation depth = %.0f \%',m(i)))

end
legend('r_i','m','cc','gzi')

SPLfm = [40 60 80];
deltaf = 700;
figure;
for i = 1:3
    yfm = fm(fc, dur, fs, fmod, deltaf); % FM_sine(fc,dur,fs,fmod,m(i),SPL);
    insigFM(:,i) = setdbspl(yfm,SPLfm(i));
    [Rfm(i) outPsy out] = Roughness_offline(insigFM(1:N,i),fs,N);
    
    subplot(3,1,i)
    plot(out.ri(1,:),'r'), hold on
    plot(out.mdept,'b-')
    plot(out.ki,'k')
    plot(out.gzi,'g--')
    grid on;
    title(sprintf('Level = %.0f',SPLfm(i)))

end
Rfmrel = Rfm / max(Rfm);
legend('r_i','m','cc','gzi')
end
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
