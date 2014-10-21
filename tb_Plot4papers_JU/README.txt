Examples in:
	Mesh.m
	Get_LPC_frames.m

N = 10;
M = 40;

    h = ImageSetup; 
    h.I_Matrix      = [N	/2,N/2];
    h.I_FontSize    = 10; 
    h.I_FontName    = 'Arial'; 
    h.I_Width       = 8;
    h.I_High        = 8;
    h.I_TitleInAxis = 1;
    h.I_Space       = [0.01,0.01];

    try
        dr = info.ylim(2)-info.ylim(1);
        if info.ylim(2) < max(max(zm))
            info.ylim(2) = max(max(zm)) + 5;
            info.ylim(1) = info.ylim(2) - dr;
        end
        h.I_Ylim = info.ylim;
    catch
        h.I_Ylim = [-105,5]; % Uncomment for fixing the limits in the y-axis
    end
    
    try
        h.I_Xlim = info.xlim;
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
    Saveas(gcf,[Get_TUe_paths('outputs') mfilename '-' info.title '-' stName]);
    close;
    close;

