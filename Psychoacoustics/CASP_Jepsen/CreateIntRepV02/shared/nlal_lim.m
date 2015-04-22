% nlal - Non-linear adaptation loops -
%        The input is divided by the output of a RC-lowpass-filter.
%        The input range is restricted to [0.00001 ... 1].
%        Impulse response invariant digital design.
%
% Usage: out = nlal_lim(in,Fs,limit,min)
%
% in         = input column vector
% fs         = sampling rate
%
% out        = output column vector
%
% See also help nlalmex.m, nlalmex.c
%
% References
%
% Dau, T. , Püschel, D. and Kohlrausch, A. (1996): "A quantitative model of the
%     `effective' signal processing in the auditory system: (I). Model structure",
%     J. Acoust. Soc. Am. 99, p. 3615-3622.
%
% Püschel, D. (1988): "Prinzipien der zeitlichan Analyse beim Hören," Doctoral Thesis,
%     Universität Göttingen
%
% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.
% $Revision: 1.00.2 beta$  $Date: 29-04-2004 16:14 $

% BM type is 1 = Gammatone, 2 = DRNL, Added by Morten 24.05.2006

% fast version SE 21.12.2010 10:14
% modified version
    % 01/06/12, C.T. Iben


function [out, divisorStage] = nlal_lim(in,Fs,limit,minValue);
% min = 2e-4;   	  %lowest signal level - initial transient solution, MSc
% SE passed as parameter
%min = 1e-7;   	  %lowest signal level - initial transient solution, MSc

len=length(in);   %signal length
out=max(in, minValue);
% SE 02/2012, included divisorStage
divisorStage = out*0;

tau1=0.005;			%first loop time-constant in sec
tau2=0.050;
tau3=0.129;
tau4=0.253;
tau5=0.500;			%last

%half
% tau1=0.002;			%first loop time-constant in sec
% tau2=0.025;
% tau3=0.65;
% tau4=0.127;
% tau5=0.250;			%last

% double
% tau1=0.01;			%first loop time-constant in sec
% tau2=0.1;
% tau3=0.260;
% tau4=0.500;
% tau5=1.0;			%last


%taul=0.02;			%overall lowpass

%--------------------------------------------------------------
% first adaption loop
b01=1/(tau1*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
a11=exp(-b01);			%a1 coefficient of the upper IIR-filter
b01=1-a11;
tmp21=sqrt(minValue);		%from steady-state relation
%---------------------------------------------------------------
% second adaption loop
b02=1/(tau2*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
a12=exp(-b02);			%a1 coefficient of the upper IIR-filter
b02=1-a12;
tmp22=minValue^(1/4);		%from steady-state relation
%---------------------------------------------------------------
% third adaption loop
b03=1/(tau3*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
a13=exp(-b03);			%a1 coefficient of the upper IIR-filter
b03=1-a13;
tmp23=minValue^(1/8);		%from steady-state relation
%---------------------------------------------------------------
% forth adaption loop
b04=1/(tau4*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
a14=exp(-b04);			%a1 coefficient of the upper IIR-filter
b04=1-a14;
tmp24=minValue^(1/16);		%from steady-state relation
%---------------------------------------------------------------
% fifth adaption loop
b05=1/(tau5*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
a15=exp(-b05);			%a1 coefficient of the upper IIR-filter
b05=1-a15;
tmp25=minValue^(1/32);		%from steady-state relation
%---------------------------------------------------------------
% overall lowpass-filter
%b0l=1/(taul*Fs);		%b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
%a1l=exp(-b0l);			%a1 coefficient of the upper IIR-filter
%b0l=1-a1l;
%tmp_l=0;			%from steady-state relation


corr = minValue^(1/32);		% to get a range from 0 to 100 model units

mult = 100/(1-corr);
%  mult = 100/(1.2*(1-corr)); % after expansion,DRNL, we need to compensate for
% "m" is added or altered by morten 26. jun 2006
if limit <=1 % m, no limitation
    
    for i=1:len
        %   tmp1=out(i);
        %if tmp1 < min
        %   tmp1=min;
        %end
        %---------------------------------------------------------------
        tmp1=out(i)/tmp21;
        tmp21 = a11*tmp21 + b01*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp22;
        tmp22 = a12*tmp22 + b02*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp23;
        tmp23 = a13*tmp23 + b03*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp24;
        tmp24 = a14*tmp24 + b04*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp25;
        divisorStage(i) = (tmp25-corr)*mult;
        tmp25 = a15*tmp25 + b05*tmp1;
        
        %--- Scale to model units ----------------------------------
        
        out(i) = (tmp1-corr)*mult;
        
        % --- LP --------
        
        %   if (lp)
        %	tmp1 = a1l*tmp_l + b0l*tmp1;
        %	tmp_l = tmp1;
        %
        %   end
        
        % out(i) = tmp1;
        
    end
else    % m, now limit
    min1 = tmp21; min2 = tmp22; min3 = tmp23;
    min4 = tmp24; min5 = tmp25;
    
    % calc values for exp fcn once
    maxvalue = (1 - min1^2) * limit - 1;
    factor1 = maxvalue * 2;
    expfac1 = -2/maxvalue;
    offset1 = maxvalue - 1;
    
    maxvalue = (1 - min2^2) * limit - 1;
    factor2 = maxvalue * 2;
    expfac2 = -2/maxvalue;
    offset2 = maxvalue - 1;
    
    maxvalue = (1 - min3^2) * limit - 1;
    factor3 = maxvalue * 2;
    expfac3 = -2/maxvalue;
    offset3 = maxvalue - 1;
    
    maxvalue = (1 - min4^2) * limit - 1;
    factor4 = maxvalue * 2;
    expfac4 = -2/maxvalue;
    offset4 = maxvalue - 1;
    
    maxvalue = (1 - min5^2) * limit - 1;
    factor5 = maxvalue * 2;
    expfac5 = -2/maxvalue;
    offset5 = maxvalue - 1;
    
    for i=1:len
        %  tmp1=out(i);
        %if tmp1 < min
        %   tmp1=min;
        %end
        %---------------------------------------------------------------
        tmp1=out(i)/tmp21;
        
        if ( tmp1 > 1 )                 % m,
            tmp1 = factor1/(1+exp(expfac1*(tmp1-1)))-offset1;  % m,
        end                             % m,
        tmp21 = a11*tmp21 + b01*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp22;
        
        if ( tmp1 > 1 )                 % m,
            tmp1 = factor2/(1+exp(expfac2*(tmp1-1)))-offset2;  % m,
        end                             % m,
        tmp22 = a12*tmp22 + b02*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp23;
        
        if ( tmp1 > 1 )                 % m,
            tmp1 = factor3/(1+exp(expfac3*(tmp1-1)))-offset3;  % m,
        end
        tmp23 = a13*tmp23 + b03*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp24;
        
        if ( tmp1 > 1 )                 % m,
            tmp1 = factor4/(1+exp(expfac4*(tmp1-1)))-offset4;  % m,
        end
        tmp24 = a14*tmp24 + b04*tmp1;
        
        %---------------------------------------------------------------
        tmp1=tmp1/tmp25;
        divisorStage(i) = (tmp25-corr)*mult;
        
        if ( tmp1 > 1 )                 % m,
            tmp1 = factor5/(1+exp(expfac5*(tmp1-1)))-offset5;  % m,
        end
        tmp25 = a15*tmp25 + b05*tmp1;
        
        %--- Scale to model units ----------------------------------
        
        out(i) = (tmp1-corr)*mult;
        
        % --- LP --------
        
        %   if (lp)
        %	tmp1 = a1l*tmp_l + b0l*tmp1;
        %	tmp_l = tmp1;
        %
        %   end
        
        %  out(i) = tmp1;
        
    end
end

% eof
