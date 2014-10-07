function out = ch_kernlaut24_two(rms, HL_ohc)
% function out = ch_kernlaut24_two(rms, HL_ohc)
%
% 1. Description:
%       Calculates main loudness for main excitation (rms) according to 
%       DIN45631 accounts for nonlinear component of hearing loss (HL_ohc)
%
%       Loudness transformation according to Chalupper2002, Equation 2
%       Data for thq and a0 from [3].
% 
% References:
%   [3] Paulus; Zwicker. "Computer Programmes for Calculating Loudness from 
%       Third-Octave Band Levels or from Critical Band Levels". ACUSTICA 
%       Vol.27 1972 issue 5 S.261 TableI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000 (new version with comments and examples on 06/01/2007)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 08/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normal hearing threshold:
thq = [42 18.5 11.5 8.3 6.7 5.5 4.8 4.3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

thq = thq+HL_ohc;

s = 10.^(0.22-0.005*[0.5:23.5])-1;
k = 0.23;
out = zeros(length(rms(:,1)),24);  

for i = 1:length(rms(:,1))
   le=rms(i,:);
   mp1= 0.04925*(1./s).^k.*10.^(0.1*k* thq);
   mp2=(1-s +s.*10.^(0.1*(le-thq))).^k -1;
   nm=mp1 .* mp2; 

   j=find(le <= thq | nm<0);   
   nm(j)=0;
   
   out(i,:)=nm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end