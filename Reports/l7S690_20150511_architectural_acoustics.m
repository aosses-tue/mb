function y = r20150511_architectural_acoustics(x)
% function y = r20150511_architectural_acoustics(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/05/2015
% Last update on: 11/05/2015 % Update this date manually
% Last use on   : 11/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 0;
Diary(mfilename,bDiary);

Lx = 4.6;
Ly = 5.8;
Lz = 3.3;
V = Lx*Ly*Lz; % 88.044;

m = [0 0 0 0 0 0 0.0022 .0066 .019]; % [m^-1]
SM = 5.76;

St12 = Lx * Lz * 2;
St34 = Ly * Lz * 2;
St56 = Lx * Ly * 2;
Stot = (St12 + St34 + St56);

T60     = [10.691	10.321	9.023	6.875	5.531	5.279	4.174	2.798	1.728]; % without
T60A    = [6.502	5.418	4.075	1.705	1.470	1.466	1.380	1.176	0.923]; % with absorption
f       = [31.5	63	125	250	500	1000	2000	4000	8000];

% 1.1 Sabine, no absorption:
alphaS1 = 0.161*V/Stot*(1./T60)-4*m*V/Stot;

% 1.2 Sabine, with absorption
Sc = Stot - SM; % S concrete
alphaS2 = 0.161*V/SM*(1./T60A)-alphaS1*Sc/SM-4*m*V/SM;

% Eyring's approach
for idx = 1:length(f);
    % tmp1    = (-.161*V - 4*m(idx)*V*T60(idx) )/(T60(idx)*Stot);
    % alpha1(idx) = 1 - exp(tmp1);
    
    tmp1 = -0.07*V/Stot/T60(idx) + 1.7391*m(idx)*V/Stot;
    alpha1(idx) = 1 - 10^(tmp1);
    
    % tmp2    = (-.161*V - 4*m(idx)*V*T60A(idx) )/(T60A(idx)*Stot);
    % alphaprom(idx) = 1 - exp(tmp2);
    
    tmp2    = -0.07*V/Stot/T60A(idx) + 1.7391*m(idx)*V/Stot;
    alphaprom(idx) = 1 - 10^(tmp2);
    
    alpha2(idx) = ( Stot*alphaprom(idx) - alpha1(idx)*Sc )/SM;
end
 
figure;
plot(alpha1,'ro--'), hold on
plot(alpha2)
xlabel('Frequency [Hz]')
ylabel('Absorption coefficient [dimenssionless]')
title('Absorption coefficient using Eyring''s formula')
h = gcf;
ha = gca;
set(ha,'XTickLabel',f);
grid on

% Saveas(gcf,'Absorption-coeff-eps','epsc');
% Saveas(gcf,'Absorption-coeff-emf','emf');

figure;
plot(alphaS1,'ro--'), hold on
plot(alphaS2)
xlabel('Frequency [Hz]')
ylabel('Absorption coefficient [dimenssionless]')
title('Absorption coefficient using Sabine''s formula')
h = gcf;
ha = gca;
set(ha,'XTickLabel',f);
grid on

figure;
plot(alpha2,'o-'), hold on
plot(alphaS2,'ro--')
xlabel('Frequency [Hz]')
ylabel('Absorption coefficient [dimenssionless]')
legend('Eyring','Sabine')
h = gcf;
ha = gca;
set(ha,'XTickLabel',f);
grid on

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
