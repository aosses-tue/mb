function [h ha] = Get_VoD_waveform_2_5(fi1,fi2,acmode,options,stPlot)
% function [h ha] = Get_VoD_waveform_2_5(fi1,fi2,acmode,options,stPlot)
%
% 1. Description:
%       
% 2. Additional info:
%   Tested cross-platform: Yes
%
% 3. Stand-alone example:
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: VoD_read_aligned
% Created on    : 31/10/2014
% Last update on: 31/10/2014 % Update this date manually
% Last use on   : 31/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

if nargin < 4
    stPlot = [];
end 

options = Ensure_field(options,'dest_folder',Get_TUe_paths('outputs'));
options = Ensure_field(options,'bSave',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
for acmode_i = acmode
    
    mode_idx = acmode_i-1;
    modus   = num2str(mode_idx); % 1 to 4
    
    % Loading audio files...
    
    [yi1, fs] = Wavread(fi1);
    [yi2, fs] = Wavread(fi2);
    
    t       = (1:length(yi1))/fs; % same t for measured and modelled

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    n = 2; 

    stPlot = Ensure_field(stPlot,'color',{'b-','r--'});
    stPlot = Ensure_field(stPlot,'LineWidth',[1 2]);

    stPlot = Ensure_field(stPlot,'Title1',sprintf('(a) ac-mode %.0f',acmode));
    stPlot = Ensure_field(stPlot,'Title2','(b) ');
    
    figure; 
    subplot(n,1,1)
    plot(t,yi1,stPlot.color{1})
    ha = gca;

    title(stPlot.Title1)
    stPlot = rmfield(stPlot,'Title1');
    ylabel('Amplitude')

    subplot(n,1,2)
    ylims = get(ha,'YLim');
    set(ha(end),'YLim',1.2*ylims); % expand YLim in 20%

    plot(t,yi2,'r-');
    ha(end+1) = gca;
    title(stPlot.Title2)
    ylims = get(ha(end),'YLim');
    set(ha,'YLim',1.2*ylims); % expand YLim in 20%

    if n == 2
        xlabel('Time [s]')
    end

    ylabel('Amplitude')

    linkaxes(ha,'x')

    h = gcf;    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end