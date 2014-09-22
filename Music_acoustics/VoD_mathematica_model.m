function VoD_mathematica_model(ac_mode,info)
% function VoD_mathematica_model(ac_mode,info)
%
% 1. Description:
%       ac_mode corresponds to the acoustic mode of resonance for the hummer.
%       Possible values are from 2 to 5.
%       PAmp = Envelope
% 
% 2. Additional info:
%   Tested cross-platform: Yes
%
% 3.1 Stand-alone example:
%       VoD_mathematica_model;
% 
% 3.2 Stand-alone example:
%       info.normalise_time_factor = 0.2688; % period, acoustic mode 5
%       VoD_mathematica_model(5,info);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 27/05/2014
% Last update on: 25/07/2014 % Update this date manually
% Last used on  : 25/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

if nargin < 2
    info = [];
end

info.bSave = 1;
info = Ensure_field(info,'normalise_time_factor',1);

paths.VoD = Get_TUe_paths('db_voice_of_dragon');

% 1. Parameters:

tanalysis = 5;
if nargin == 0
    ac_mode = 2; % from 2 to 5
end

fieldtype = 2; % 1 = far-field, 2 = near-field

% Trot = 344/1000;
params      = Get_VoD_params(0);

if fieldtype == 2 % Near field
    Xo  = params.Xo_nf
    Yo  = params.Yo_nf;
    Zo  = params.Zo_nf;
else
    Xo  = params.Xo_ff;
    Yo  = params.Yo_ff;
    Zo  = params.Zo_ff;
end
Pos_obs = [Xo, Yo, Zo]; % Observer's position (microphone)

mode_idx    = ac_mode - 1;
Trot        = params.Tmodel(mode_idx);
mf          = params.mf(mode_idx);
fn          = params.mf;
% Trotation   = params.Tmodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nakiboglu2012:
n   = [2 5]; % modes

cf  = 1.78e-2; % resistant coeff, experimentally determined
co = 340; % speed of sound [m/s]
ratio_Vin_Vtot = 0.83; % page 753
W   = 5e-3; % see Nakiboglu2011
rup = 1e-3; % see Nakiboglu2011
Srwrup = 0.44; % Strouhal number, page 753
L       = 70/100; % check if 75 cm is better
Din     = 26.5e-3; % Inner diameter [m], page 751
ceff    = co*sqrt(ratio_Vin_Vtot); % effective speed of sound [m/s], approx 310

% Leff = mean(n*ceff./(2*fn)); % approx. 0.7264 m, using Eq.10, pag 753

Omegan = fn*(W+rup)/(Srwrup*params.R)*sqrt(1+4*cf*L/Din);
2*pi./Omegan % Periods

% fn.*Trotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lp = 70/100;
Sp = pi/4*(33/100)^2;
Rp = params.R;

n_exact = 2*Lp*mf/ceff;
n   = round(n_exact); 

wn  = pi*n_exact*ceff/Lp; % Mathematica model: wn
%   - wn must be related with modal frequeny as: wn = 2*pi*mf
%   - To fulfill this, n_exact has to be used. It is an additional step though,
%     since mf is known...

Omegao = 2*pi/Trot;
fs  = 44100;
t   = 0:1/fs:tanalysis; %Trot; % 0 to 1 rotation period
Z1 = 220/100;

S1r = [0 0  Z1]; % Real position, fixed source
S1m = [0 0 -Z1]; % Im   position, fixed source
S2r = [Rp*cos(Omegao*t); Rp*sin(Omegao*t);   Z1+0*t]; % Real position, rotating source
S2m = [Rp*cos(Omegao*t); Rp*sin(Omegao*t);  -Z1+0*t]; % Im   position, rotating source

figure;
plot3(S2r(1,1), S2r(2,1),S2r(3,1),'x','LineWidth',12), hold on
plot3(Xo, Yo, Zo, 'rx','LineWidth',12) % Mic position
plot3(S2r(1,:),S2r(2,:),S2r(3,:))
grid on
xlabel('X position [m]')
ylabel('Y position [m]')
legend('Starting VoD point','Mic position')

Pos_obsg = repmat(Pos_obs',1,size(S2r,2));
Dist = sqrt( sum((S2r - Pos_obsg).^2) );
figure;
plot(t/info.normalise_time_factor,Dist)
ylabel('Distance rotating source to mic. [m]')
xlabel('Time [s]')
grid on

rhoo = 1;

% 3. Sound radiation without doppler:

Qn  = Sp/(ceff*rhoo) * cos(wn*t); % Mathematica model: Q1

X   = [Xo*cos(Omegao*t); Xo*sin(Omegao*t);   Z1+0*t];
X1  = repmat([0 0 Z1]',1,length(X));
X2  = repmat([Rp 0 Z1]',1,length(X));
XminX1 = sqrt(sum((X-X1).*(X-X1)));
XminX2 = sqrt(sum((X-X2).*(X-X2)));

p   = rhoo*wn*Qn/(4*pi);
Amplitude = abs( exp(-j*wn*(XminX1)/co)./XminX1 + (-1)^n * exp(-j*wn*(XminX2)/co)./XminX2 );
pTotal_no_doppler = p.*Amplitude;

figure;
plot(t,pTotal_no_doppler,'b',t,Amplitude,'r');

figure;
freqz(pTotal_no_doppler,1,4096*2)

filename = sprintf('%smodel-ac-mode-%.0f-no-doppler',Get_TUe_paths('outputs'),n);
Wavwrite(pTotal_no_doppler,fs,filename)

% 4. Sound radiation with doppler:

te  = t - XminX2/co;
DSr2 = [-Rp*Omegao*sin(Omegao*t);Rp*Omegao*cos(Omegao*t);0*t];
Ms2 = (X-X2)./ repmat(XminX2,3,1) .* DSr2/co;
% pd  = -rhoo*wn*Qn/(4*pi)*( sin(t-XminX1/co)./XminX1 + (-1)^(n+1)*sin(wn*te))/(1-Ms2) ;
% 
% figure;
% plot(t,pd )% ,'b',t,Amplitude,'r');
% 
% Wavwrite(pd,fs,[Get_TUe_paths('outputs') 'test-doppler'])

% 
% ti = 0;
% % tf = 2* (2*pi)/w;
% 
% t = ti: 1/srate : tfin - 1/srate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end