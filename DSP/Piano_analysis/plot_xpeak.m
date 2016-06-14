function pipo = plot_xpeak(APEAK,FPEAK,FE,SEUIL)

%usage
%    function plot_xpeak(APEAK,FPEAK,FE,SEUIL)
% les amplitudes APEAK sont en lineaire

%AXE DES X
FF = [FPEAK;FPEAK];
%AXE DES Y
MAX=max(APEAK);
AA = [ones(size(1:length(APEAK)))*(20*log10(MAX+eps) + SEUIL) ;20*log10(APEAK+eps)];


MAX=ceil(max(20*log10(APEAK+eps))/10)*10;
MIN = MAX + SEUIL;


plot(FF/1000,AA,'-'); drawnow;
F_max=FE/2000;
axis([ 0 F_max MIN MAX]);
%axis([ 0 20 MIN MAX]);
title(' Detection de Pics');
if FE==1
xlabel('Frequence reduite')
else
xlabel('Frequence (kHz)')
end
ylabel('Amplitude (dB)')
axis;
