% mfbtd.m - applies modulation filterbank	as suggested in Dau et al. 1997.
%
% Usage: [inf, out] = mfbtd(in,lmf,umf,style,fs)
%
% in    = input column vector.
% lmf   = lowest modulation-filter center frequency,
%         if 0 the output of a 2nd-order
%         Butterworth lowpass filter is additionally computed.
% umf   = highest modulation-filter center frequency,
%         for typical applications choose umf = 1500.
%         If lmf = umf only the output of a single filter
%       	 at lmf is computed.
% style = determines fc of the lowpass filter: 1 = 2.5 Hz, 2 = 7.5 Hz.
% fs    = sampling rate in Hz,
% 		    should be greater than 8000 Hz to avoid aliasing errors.
%
% inf							= center frequencies of the modulation filters.
% [out1,out2, ...,outn] = each column of martrix out contains the output of
%								  a single modulation filter.
%
% copyright (c) 1999 Stephan Ewert and Torsten Dau, Universitaet Oldenburg

function [inf, out] = mfbtd(in,lmf,umf,style,fs)

if fs < 8000
    warning('sample rate lower than 8000 Hz')
end

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));
%lpcut = 150;
%ex=((3+sqrt(5))/2)^(1/den);

if lmf == 0
    sw = 1;
    switch style
        case 1
            startmf = 5;
            b2lpcut = 2.5;
            wb2lp = 2*pi*b2lpcut/fs;
            [b2,a2] = solp(wb2lp,1/sqrt(2));
        case 2
            startmf = 10;
            b2lpcut = 7.5;
            wb2lp = 2*pi*b2lpcut/fs;
            [b2,a2] = solp(wb2lp,1/sqrt(2));
    end
elseif lmf > 0
    startmf = lmf;
    sw = 2;
end

if lmf == umf
    mf = startmf;
    sw = 0;
    if lmf == 0
        sw =3;
    end
else
    if startmf >= 10
        tmp = fix(log(umf/startmf)/log(ex));
        mf = 0:tmp;
        mf = ex.^mf;
        mf = startmf*mf;
    else
        tmp = fix((min(umf,10) - startmf)/bw); %changed
        tmp = 0:tmp;
        mf = startmf + 5*tmp;
        tmp2 = (mf(end)+bw/2)/(1-1/(2*Q));
        tmp = fix(log(umf/tmp2)/log(ex));
        tmp = 0:tmp;
        tmp = ex.^tmp;
        mf=[mf tmp2*tmp];
    end
end


%[b2,a2] = folp(wb2lp);

switch sw
    case 0									% only one modulation filter
        w0 = 2*pi*mf/fs;
        if mf < 10
            [b3,a3] = efilt(w0,2*pi*bw/fs);
        else
            [b3,a3] = efilt(w0,w0/Q);
        end
        out = 2*filter(b3,a3,in);
        inf = mf;
    case 1									% lowpass and modulation filter(s)
        out = zeros(length(in),length(mf)+1);
        out(:,1) = filter(b2,a2,in);
        for i=1:length(mf)
            w0 = 2*pi*mf(i)/fs;
            if mf(i) < 10
                [b3,a3] = efilt(w0,2*pi*bw/fs);
            else
                [b3,a3] = efilt(w0,w0/Q);
            end
            out(:,i+1) = 2*filter(b3,a3,in);
        end
        inf = [0 mf];
    case 2									% only modulation filters
        out = zeros(length(in),length(mf));
        for i=1:length(mf)
            w0 = 2*pi*mf(i)/fs;
            if mf(i) < 10
                [b3,a3] = efilt(w0,2*pi*bw/fs);
            else
                [b3,a3] = efilt(w0,w0/Q);
            end
            out(:,i) = 2*filter(b3,a3,in);
        end
        inf = mf;
    case 3									% only lowpass
        out = filter(b2,a2,in);
        inf = 0;
end

% subfunctions

% complex frequency shifted first order lowpass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];

% second order Butterworth lowpass filter
function [b,a] = solp(w0,Q)

W0 = tan(w0/2);

b = [1; 2; 1];
a = [1 + 1/(Q*W0) + 1/W0^2; -2/W0^2 + 2; 1 - 1/(Q*W0) + 1/W0^2];

b = b/a(1);
a = a/a(1);

% eof






