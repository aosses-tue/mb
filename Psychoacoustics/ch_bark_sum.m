function kl= ch_bark_sum(ns)
% function kl= ch_bark_sum(ns)
%
% spectral summation of specific loudness in 0.1 Bark steps to 1 Bark steps
% kl is a Matrix with Dimension [t(1),t(2)/10] ; t=size(ns)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000 (new version with comments and examples on 06/01/2007)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 08/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=size(ns);

for i=1:t(1) 
   for j=1:t(2)/10 
      kl(i,j)=0.1*sum(ns(i,(10*j-9):(10*j)));
   end
end

ind=find(kl==0);
kl(ind)=realmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
