function h = Mesh(t,f,z,options)
% function h = Mesh(t,f,z,options)
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
% 2. Stand-alone example:
%       x = 0:.2:8;
%       y = 0:10:5000;
%       z = rand(length(y),length(x));
%       Mesh(x,y,z);
% 
% 3. Additional info:
%   Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 23/07/2014
% Last update on: 21/08/2014 % Update this date manually
% Last used on  : 02/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    options = [];
end

options = Ensure_field(options,'bPlot3D',1);
options = Ensure_field(options,'bPlot2D',~options.bPlot3D);
options = Ensure_field(options,'title','');
options = Ensure_field(options,'step1',1);
options = Ensure_field(options,'step2',1);

step_t = options.step1;
step_f = options.step2; % round(length(f)/M);

N = round(length(t)/step_t);
M = round(length(f)/step_f);

[tm fm] = meshgrid(t(1:step_t:end),f(1:step_f:end));

zm = z(1:step_f:end, :);
zm = zm(:,1:step_t:end);

if options.bPlot3D
    meshz(tm,fm,zm)
    % plot3(tm,fm,zm)
    
    xlabel('time')
    ylabel('frequency')
    zlabel('F(t,f)')
end
 
if options.bPlot2D
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
        
        if i == 1
            try
                title(options.title);
            end
        end
    end

    h = ImageSetup; 
    h.I_Matrix      = [N/2,N/2];
    h.I_FontSize    = 10; 
    h.I_FontName    = 'Arial'; 
    h.I_Width       = 8;
    h.I_Height      = 8;
    h.I_TitleInAxis = 1;
    h.I_Space       = [0.01,0.01];

    try
        dr = options.ylim(2)-options.ylim(1);
        if options.ylim(2) < max(max(zm))
            options.ylim(2) = max(max(zm)) + 5;
            options.ylim(1) = options.ylim(2) - dr;
        end
        h.I_Ylim = options.ylim;
    catch
        h.I_Ylim = [-105,5]; % Uncomment for fixing the limits in the y-axis
    end
    
    try
        h.I_Xlim = options.xlim;
    end
    
    % h.I_Xlim = [0,5];
    h.I_Grid = 'on'; 
    h.I_KeepColor = 0; 
    h.I_Handles = gcf;
    h.prepareAllFigures;
    h.arrayAddedHandles = 1;
    add2ArraySubplotVer(h);

    stName = Get_date;
    stName = stName.date2print;
    Saveas(gcf,[Get_TUe_paths('outputs') mfilename '-' options.title '-' stName]);
    close;
    close;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
