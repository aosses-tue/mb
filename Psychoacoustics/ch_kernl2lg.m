function lg = ch_kernl2lg(kernl,HL_ohc,HL_ihc,fg)
% function lg = ch_kernl2lg(kernl,HL_ohc,HL_ihc,fg)
% 
% 1. Description:
%   Transforms critical band loudness into critical band level inverse 
%   loudness transformation assumes 24x1 vectors if less than 4 arguments  
%
% required functions: korrel.m, kernl2lg.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000
% Edited on     : 06/02/2007 (new version with comments and examples)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 29/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin==1)
   HL_ohc=zeros(1,24);
   HL_ihc=zeros(1,24);
end   


s=10.^(0.22-0.005*[0.5:23.5])-1;
k=0.23;

a0 = [ 0 0 0 0 0 0 0 0 0 0 -.2 -.5 -1.2 -2.1 -3.2 -4.6 -5.5 -5.6 -4.3 -2.5 -0.1 2.8 6.4 20.0];

% Normal hearing threshold:
thq=[42 18.5 11.5 8.3 6.7 5.5 4.8 4.3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

if nargin==4
   a0=a0(fg);
   thq=thq(fg);
end   

lhs=thq+HL_ohc; 

c1=0.0533*(1./s).^k.*10.^(0.1.*k.*lhs);
exact=((kernl./c1+1).^(1/k) -1 +s)./s;
lg=10.*log10(exact)+lhs;

lg=lg+a0+HL_ihc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end