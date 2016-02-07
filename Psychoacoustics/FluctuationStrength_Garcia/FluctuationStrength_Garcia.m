function [fluct fi] = FluctuationStrength_Garcia(insig_b, fs, N, model_par)
% function [fluct fi] = FluctuationStrength_Garcia(insig_b, fs, N, model_par)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
% 
%       Comments on this implementation: cross-correlation seems not to be 
%       working properly.
%
%       Adapted from Model/Helper/FluctuationStrength.m (first adaptation from Roughness model)
%       Adapted from Model/Helper/FluctuationStrength_debug.m (main file as in the thesis)
% 
%       Comments:
%           fs  - should be an input parameter
%           N   - should be an input parameter
%           Too difficult to change N!
%           gzi - where was it taken from?
% 
% 2. Stand-alone example:
%        
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Rodrigo Garcia L./Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Previous name : FluctuationStrength_Garcia_offline.m
% Created on    : 03/02/2016
% Last update on: 03/02/2016 
% Last use on   : 03/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    model_par = Get_fluctuation_strength_params(N,fs);
    model_par.debug = 'none';
end

model_par = ef(model_par,'window_type','blackman');

% FluctuationStrength_debug.m

%% ei = peripheral_stage(insig,fs,N);
% Blackman window
switch model_par.window_type
    case 'blackman'
        window = blackman(N);
        window = window/mean(window);
    case 'cosine'
        window = ones(N,1);
        attackrelease = 50;
        window = Do_cos_ramp(window,fs,attackrelease);
        % window = From_dB(100)*window; %/mean(window); % No correction to the window
end

overlap = round(0.9*N);
insig_b = buffer(insig_b,N,overlap,'nodelay');
nFrames = size(insig_b,2);

fluct   = zeros(1,nFrames); % Memory allocation

for iFrame = 1:nFrames
    
    insig = insig_b(:,iFrame);
    % Apply window to frame
    insig = transpose(window .* insig);

    %% 1. Peripheral stages

    % 1.1 Peripheral hearing system
    insig = PeripheralHearingSystem(insig,fs,N);

    % 1.2 Excitation patterns
    dBFS = 100; % unit amplitude corresponds to 100 dB (AMT Toolbox convention)
    ei = TerhardtExcitationPatterns(insig,fs,dBFS);
    %%%

    %% 2. Modulation depth (estimation)
    [mdept,hBPi,ei] = il_modulation_depths(ei,model_par.Hweight);
    %%%    

    %% 3. Cross-correlation coefficient:
    switch model_par.dataset
        case {0,99}
            Ki = il_cross_correlation(abs(ei));
            fi = il_specific_fluctuation(mdept,Ki,model_par,model_par.dataset);

        case 1
            Ki = il_cross_correlation(hBPi);
            fi = il_specific_fluctuation(mdept,Ki,model_par,model_par.dataset);
    end

    fi        = model_par.cal * fi;
    fluct(iFrame) = sum(fi);
end

disp('')
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdept,hBPi,ei] = il_modulation_depths(ei,Hweight)
    [Chno,Nc] = size(ei);

    hBPi   = zeros(Chno,Nc);
    hBPrms = zeros(1,Chno);
    mdept  = zeros(1,Chno);

    ei      = abs(ei);
    h0      = mean(ei,2);
    ei      = ei - repmat(h0,1,Nc);
    for k = 1:Chno
        if ~isnumeric( Hweight )
            hBPi(k,:) = filter(Hweight,ei(k,:));
        else
            hBPi(k,:) = sosfilt(Hweight,ei(k,:));
        end
    
        hBPrms(k) = rms(hBPi(k,:));

        if h0(k) > 0
            mdept(k) = hBPrms(k)/h0(k);
            if mdept(k) > 1
                mdept(k) = 1;
            end
        else
            mdept(k) = 0;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ki = il_cross_correlation(hBPi)
    [chno,~] = size(hBPi);

    ki = zeros(1,chno-2);
    for k=1:chno-2
        cfac = cov(hBPi(k,:),hBPi(k+2,:));
        den  = diag(cfac);
        den  = sqrt(den*den');

        if den(2,1) > 0
            ki(k) = cfac(2,1)/den(2,1);
        else
            ki(k) = 0;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fi = il_specific_fluctuation(mdept,Ki,model_par,dataset)
    
    gzi = model_par.gzi; 
    p_g = model_par.p_g;
    p_m = model_par.p_m;
    p_k = model_par.p_k;
    
    Chno = length(gzi);

    fi = zeros(1,Chno);
    for k = 1:Chno
        
        md(k) = mdept(k)-0.1;
        if md(k) < 0
            md(k) = 0;
        end

        if k == Chno-1 || k == Chno
            ki = 1;
        else
            ki = Ki(k);
        end

        if k == 1 || k == 2
            ki2 = 1;
        else
            ki2 = Ki(k-2);
        end
        kp(k) = abs(ki*ki2);

    end
    
    switch dataset
        case {0,99}
            fi = gzi.^p_g .* mdept.^p_m .* kp.^p_k;
            pgtest = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 3 4];
            calval = model_par.cal;
            for i = 1:length(pgtest)
                fitest(i,1:47) = gzi.^pgtest(i) .* mdept.^p_m .* kp.^p_k;
                FStest(i) = calval*sum(fitest(i,:));
            end
        case 1
            fi = gzi.^p_g*md.^p_m*kp.^p_k;
    end

    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%     for k = 1:Chno
%         
%         % Fei         = fft( etmp_td(k,:) - h0(k) );
%         % hBPi(k,:)   = 2 * real(ifft(Fei.* Hweight'));
%         hBPi(k,:) = 2*filter(Hweight,etmp_td(k,:)-h0(k));
%         
%     end
%     
%     hBPrms  = rms( hBPi ,'dim',2);
%     
%     idx = find(h0 > 0);
%     mdept(idx) = hBPrms(idx) ./ h0(idx);
%     
%     idx = find(mdept > 1);
%     mdept(idx) = 1;
%     
%     %%%
%     warning('Temporal')
%     Md = mdept - 0.1;
%     idx = find(Md < 0);
%     Md(idx) = 0;
%     %%%
%     
%     ki = zeros(1,Chno - 2);
%     fi = zeros(1,Chno);
% 
%     % Find cross-correlation coefficients
%     for k=1:1:Chno-2
%         cfac    = cov(hBPi(k,:),hBPi(k + 2,:));
%         den     = diag(cfac);
%         den     = Round( sqrt(den * den'), 10); % rounded to 10 decimals
% 
%         if den(2,1) > 0
%             ki(k) = cfac(2,1) / den(2,1);
%         elseif den(2,1) < 0
%             ki(k) = 0;
%         else
%             ki(k) = 0;
%         end
%     end
% 
%     % Calculate specific fluctuation strength fi and total FS
%     fi(iFrame,1) = gzi(1)^p_g * Md(1)^p_m * abs(ki(1))^p_k;
%     fi(iFrame,2) = gzi(2)^p_g * Md(2)^p_m * abs(ki(2))^p_k;
% 
%     for k = 3:1:45
%         fi(iFrame,k) = gzi(k)^p_g * Md(k)^p_m * abs( ki(k - 2) * ki(k) )^p_k;
%     end
% 
%     fi(iFrame,46) = gzi(46)^p_g * Md(46)^p_m * abs(ki(44))^p_k;
%     fi(iFrame,47) = gzi(47)^p_g * Md(47)^p_m * abs(ki(45))^p_k;
% 
%     FS(iFrame) = Cal * sum(fi(iFrame,:));
% end
% 
% dataOut{1} = FS;
% dataOut{2} = fi;
% dataOut{3} = SPL;
% 
% nParam      = 1;
% out.Data1   = transpose(FS);
% output.name{nParam} = 'Fluctuation strength';
% output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end
% 
% function gzi = il_create_gzi(Chno);
% 
%     Chstep = 0.5;
%     
%     g0 = [
%         0       0
%         125     1
%         250     1
%         500     1
%         1000    1
%         1500    1
%         2000    1
%         3000    1
%         4000    1
%         6000    1
%         8000    1
%         16000   0
%     ];
%     
%     gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
%     gzi(isnan(gzi)) = 0;
%     
%     gzi = ones(1,Chno);