function [outsig, fc] = drnl_CASP(insig,fs,BM)
% function [outsig, fc] = drnl_CASP(insig,fs,BM)
% 
% script for pemo preprocessing using an implementation of the dual
% resonance nonlinear (DRNL) filter (Lopez-Poveda, meddis 2001)
% The filter models the BM non-linearity
% Author: Morten Loeve Jepsen, 2.nov 2005
%
% usage: out = drnl(x,CF,fs)
% Original file: '..\Psychoacoustics\CASP_Jepsen\CreateIntRepV02\casp2008\bm\drnl.m'
% Created on  : 24/04/2015
% Last used on: 13/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    % necessary parameters for basilar-membrane and modulation filterbank
    [BM MF]     = CaspPreProcCfg;
    [BM MF Lp]  = CaspPreProcInit(BM, MF, fs);
    warning('In this implementation IntRep neither MF are not being used...')
end

fc = BM.CenterFreqs;
NrChannels = length(fc);
outsig = zeros(length(insig),NrChannels);

for ii = 1:NrChannels
    
    CF = fc(ii); % current centre frequency
    [linDRNLpar,nlinDRNLpar] = getDRNLparam(CF);
                                 %(fc                , BW               ,fs) 
    [GTlin_b,GTlin_a] = coefGtDRNL(linDRNLpar(5).vals,linDRNLpar(3).vals,fs); %get GT filter coeffs
                                 %(fc                ,fs)
    [LPlin_b,LPlin_a] = coefLPDRNL(linDRNLpar(5).vals,fs); % get LP filter coeffs

    y_lin = insig.*linDRNLpar(4).vals; % Apply linear gain

    % Now filtering
    for n = 1:linDRNLpar(2).vals % Gammatone filtering multiple times for cascading
        y_lin = real(filter(GTlin_b,GTlin_a,y_lin));
    end
    for n = 1:linDRNLpar(6).vals % cascade of lowpass filters
        y_lin = filter(LPlin_b,LPlin_a,y_lin);
    end
    % end of linear part %%%%%%%%%%%%%%%%%%%%%%%

    % Non-linear part%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [GTnlin_b,GTnlin_a] = coefGtDRNL(nlinDRNLpar(7).vals,nlinDRNLpar(3).vals,fs); %get GT filter coeffs
    [LPnlin_b,LPnlin_a] = coefLPDRNL(nlinDRNLpar(7).vals,fs); % get LP filter coeffs

    y_nlin = insig;

    % Now GT filtering
    for n = 1:nlinDRNLpar(2).vals % Gammatone filtering multiple times for cascading
        y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));
    end

    % Broken stick nonlinearity
    a = nlinDRNLpar(4).vals;
    b = nlinDRNLpar(5).vals;
    c = nlinDRNLpar(6).vals;

    y_decide = transpose( [a*abs(y_nlin) b*(abs(y_nlin)).^c] ); 

    % SE 10.05.2010 12:39 WTF?
    % pretty costy
    %[tmp,tmpindex] = min(y_decide);
    %fid = fopen('display.txt','a');
    %fwrite(fid,num2str(tmpindex));
    %fclose(fid);

    y_nlin = transpose( transpose( sign(y_nlin) ).* min(y_decide) );

    % Now GT filtering again
    for n = 1:nlinDRNLpar(2).vals % Gammatone filtering multiple times for cascading
        y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));
    end
    % then LP filtering
    if nlinDRNLpar(8).vals > 0
        for n = 1:nlinDRNLpar(8).vals % cascade of lowpass filters
            y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);
        end
    end

    outsig(:,ii) = (y_lin + y_nlin);
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eof

