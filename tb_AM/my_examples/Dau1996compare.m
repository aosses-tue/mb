function outs = Dau1996compare(insig1,insig2,fs,opts)
% function outs = Dau1996compare(insig1,insig2,fs,opts)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 0;
%       demo_dau1996b(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 08/10/2014
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 22/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

t = ( 1:length(insig1) )/fs;
t = t(:);

if nargin < 4
    opts = [];
    if nargout == 0
        opts.bPlot = 1;
    else
        opts.bPlot = 0;
    end
end

opts = Ensure_field(opts,'fc_idx',3000);
bPlot = opts.bPlot;
 

%% Processing 

insig2 = insig1 + insig2;

[out1, fc , outsig1] = dau1996preproc(insig1,fs);
[out2, fc2, outsig2] = dau1996preproc(insig2,fs);

idx     = max(find(fc<opts.fc_idx));
fcentre = fc(idx);

haxis = [];

% Normalisation of the template:
% outs.template = outsig2.out04_LPF(:,idx)-outsig1.out04_LPF(:,idx);
outs.template_no_norm   = outsig2.out04_LPF-outsig1.out04_LPF;

outs.template       = Normalise_signal(outs.template_no_norm);

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
end

outs.t = t;
outs.h = h;
outs.outsig1 = outsig1;
outs.outsig2 = outsig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
