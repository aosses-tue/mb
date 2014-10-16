function outs = Dau1996compare(insig1,insig2,fs,bPlot)
% function outs = Dau1996compare(insig1,insig2,fs,bPlot)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 08/10/2014
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 15/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

t = ( 1:length(insig1) )/fs;
t = t(:);

if nargin < 4
    if nargout == 0
        bPlot = 1;
    else
        bPlot = 0;
    end
end

%% Processing 

insig2 = insig1 + insig2;

[out1, fc ,outsig1] = dau1996preproc(insig1,fs);
[out2, fc2,outsig2] = dau1996preproc(insig2,fs);

idx     = max(find(fc<3000));
fcentre = fc(idx);

haxis = [];

% outs.template = outsig2.out04_LPF(:,idx)-outsig1.out04_LPF(:,idx);
outs.template_no_norm   = outsig2.out04_LPF-outsig1.out04_LPF;
outs.template           = outs.template_no_norm ./ repmat( abs(sum(outs.template_no_norm)), size(outs.template_no_norm,1),1) ;
outs.idx                = idx; % plotted/to be plotted idx

if bPlot
    
    K = 4096*2;
    options.fs = fs;
    
    %% Similar to Dau1996a, Figure 6
    % Apparently this filter bank is the replaced one (see Dau1997b)
    figure;
    freqfft(outsig1.out01_filterbank(:,idx) ,K,options);
    % legend(title1)
    xlim([0 8000])
    ylim([-30 20])
    xlabel('Time [s]')
    ylabel('Amplitude')
    h(end+1)=gcf;

    %% Similar to Dau1996a, Figure 7
    for i = idx 

        figure
        subplot(3,1,1)
        plot(   t, outsig1.out04_LPF(:,idx) ); hold on

        subplot(3,1,2)
        plot(   t, outsig2.out04_LPF(:,idx) ); hold on

        subplot(3,1,3)
        plot(   t, outs.template(:,idx) ); hold on

        h(end+1)=gcf;
        haxis = gca;
    end
    linkaxes(haxis,'x');
    % xlim([0 0.5])

    % figure; 
    % mesh(outsig1.out03_adaptloop)
    % xlabel('Band number'); 
    % ylabel('Sample number'); 
    % zlabel('Amplitude [MU]'); 
    % h(end+1)=gcf;
    % ha = gca;

    % zlim([yminMU ymaxMU])
    % set(ha,'CameraPosition',[79.9938 138832 22702.7]);
    % colorbar('vert')
    % 	CameraPosition = [79.9938 138832 22702.7]
end

outs.t = t;
outs.h = h;
outs.outsig1 = outsig1;
outs.outsig2 = outsig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
