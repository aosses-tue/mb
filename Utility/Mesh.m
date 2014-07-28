function h = Mesh(t,f,z,info)
% function h = Mesh(t,f,z,info)
%
% 1. Description:
%       It generates a 3D plot using a reduced time and frequency resolution. 
%       By default the time axis is sampled at N = 10 times from tmin to tmax.
%       By default the frequency axis is sampled at M = 40 times from fmin to fmax
%       In addition, a figure consisting of 10 2D-subplots is generated and
%       stored, where each subplot corresponds to ti and level is plotted as
%       function of frequency. This plot is automatically closed, but it's
%       being stored at 'outputs' directory.
% 
% 2. Additional info:
%   Tested cross-platform: YES
%  
% 3. Stand-alone example:
%       x = 0:.2:8;
%       y = 0:10:5000;
%       z = rand(length(y),length(x));
%       Mesh(x,y,z);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/07/2014
% Last update on: 28/07/2014 % Update this date manually
% Last used on  : 28/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    info = [];
end

info = Ensure_field(info,'bPlot3D',1);
info = Ensure_field(info,'bPlot2D',~info.bPlot3D);

N = 10;
M = 40;
step_t = round(length(t)/N);
step_f = round(length(f)/M);

[tm fm] = meshgrid(t(1:step_t:end),f(1:step_f:end));

zm = z(1:step_f:end, :);
zm = zm(:,1:step_t:end);

if info.bPlot3D
    mesh(tm,fm,zm)
    xlabel('time')
    ylabel('frequency')
    zlabel('F(t,f)')
end
 
if info.bPlot2D
    figure;
    tm = t(1:step_t:end);
    fm = f;

    [tm fm] = meshgrid(t(1:step_t:end),f);
    zm = z(:,1:step_t:end);

    for i = 1:N
        subplot(N/2,2,i);
        plot(fm(:,i),zm(:,i))
        title(sprintf('t = %.2f [s]',tm(1,i)))
        xlabel('Freq. [Hz]')
        ylabel(sprintf('Amplitude\n[dBFS]'))
    end

    h = ImageSetup; 
    h.I_Matrix      = [N/2,N/2];
    h.I_FontSize    = 10; 
    h.I_FontName    = 'Arial'; 
    h.I_Width       = 8;
    h.I_High        = 8;
    h.I_TitleInAxis = 1;
    h.I_Space       = [0.01,0.01];

    h.I_Ylim = [-105,5]; % Uncomment for fixing the limits in the y-axis
    % h.I_Xlim = [0,5];
    h.I_Grid = 'on'; 
    h.I_KeepColor = 0; 
    h.I_Handles = gcf;
    h.prepareAllFigures;
    h.arrayAddedHandles = 1;
    add2ArraySubplotVer(h);

    stName = Get_date;
    stName = stName.date2print;
    Saveas(gcf,[Get_TUe_paths('outputs') mfilename '-' stName]);
    close;
    close;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end