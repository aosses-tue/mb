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
options = Ensure_field(options,'Title','');
options = Ensure_field(options,'step1',1);
options = Ensure_field(options,'step2',1);
options = Ensure_field(options,'ZLabel','Amplitude');
options = ef(options,'XLabel','x-axis');
options = ef(options,'XLabel','y-axis');

step_t = options.step1;
step_f = options.step2; % round(length(f)/M);

bPlot2D = options.bPlot2D;
bPlot3D = options.bPlot3D;
ZLabel  = options.ZLabel;
[tm fm] = meshgrid(t(1:step_t:end),f(1:step_f:end));

zm = z(1:step_f:end, :);
zm = zm(:,1:step_t:end);

if bPlot3D
    meshz(tm,fm,zm)
    % meshz(tm,fm,zm,gradient(zm))
    % plot3(tm,fm,zm)
    % colormap([0.9 0.9 0.9])
    xlabel(options.XLabel)
    ylabel(options.YLabel)
    zlabel(options.ZLabel)
end
 
if bPlot2D
    
    Ntmp = round(length(t)/step_t);
    Mtmp = round(length(f)/step_f);

    N = min(Ntmp,Mtmp);
    M = max(Ntmp,Mtmp);

    if bPlot2D == 1
        options = Ensure_field(options,'I_Matrix',[N/2 N/2]);
        options = Ensure_field(options,'I_Width' ,8);
        options = Ensure_field(options,'I_Height',8);
        
        if options.I_Matrix(1)*options.I_Matrix(2) < N
            warning('Resizing I_Matrix')
            options.I_Matrix = [N/2 N/2];
        end
    end
    
    figure;
    tm = t(1:step_t:end);
    fm = f;

    [tm fm] = meshgrid(t(1:step_t:end),f);
    zm = z(:,1:step_t:end);

    if Ntmp>Mtmp
        bPlot_time = 1;
        bPlot_freq = 0;
        fm = transpose(fm);
        zm = transpose(zm);
        tm = transpose(tm);
        XLabel = options.XLabel;
    else
        bPlot_time = 0;
        bPlot_freq = 1;
        XLabel = options.YLabel;
    end
       
    ha = [];
    MinAmp = 0;
    MaxAmp = 0;
    
    for i = 1:N
        subplot(N/2,2,i);
        z2plot = zm(:,i);
        MinAmp = min(MinAmp,min(z2plot));
        MaxAmp = max(MaxAmp,max(z2plot));
        if bPlot_freq
            plot(fm(:,i),z2plot)
            strtitle = sprintf('t = %.2f [s]',tm(1,i));
            title(strtitle);
            ha(end+1) = gca;
        elseif bPlot_time
            plot(tm(:,i),z2plot)
            strtitle = sprintf('f = %.1f [Hz]',fm(1,i));
            title(strtitle);
            ha(end+1) = gca;
        end
        ylabel(ZLabel)
        xlabel(XLabel)
        
        if i == 1
            try
                title(sprintf('%s\n%s',options.Title,strtitle));
            end
        end
    end
    linkaxes(ha,'xy');
    ylim([MinAmp MaxAmp]);
    
    h = ImageSetup; 
    h.I_Matrix      = options.I_Matrix;
    h.I_FontSize    = 10; 
    h.I_FontName    = 'Arial'; 
    h.I_Width       = options.I_Width;
    h.I_Height      = options.I_Height;
    h.I_TitleInAxis = 1;
    h.I_Space       = [0.01,0.01];

    try
        h.I_Ylim = options.I_YLim;
    end
    
    h.I_Grid = 'on'; 
    h.I_KeepColor = 0; 
    h.I_Handles = gcf;
    h.prepareAllFigures;
    h.arrayAddedHandles = 1;
    
    if bPlot_freq
        add2ArraySubplotVer(h);
    elseif bPlot_time
        add2ArraySubplotHor(h);
    end

    stName = Get_date;
    stName = stName.date2print;
    Saveas(gcf,[Get_TUe_paths('outputs') mfilename '-' name2figname(options.Title) '-' stName]);
%     close;
%     close;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
