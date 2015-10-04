function [y noise] = Add_gaussian_noise(x,mu,sigma)
% function [y noise] = Add_gaussian_noise(x,mu,sigma)
%
% 1. Description:
%       It adds gaussian noise (mean mu, std deviation sigma) to the input 
%       signal x. The array x should be a column vector.
%       The variance of the noise is var = sigma.^2;
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 04/08/2015
% Last update on: 02/10/2015 
% Last use on   : 02/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1);
M = size(x,2);

noise = normrnd(mu,sigma,N,M);
y = x + noise;

if nargout == 0
    [pdfi xi] = Add_pdf2plot(noise);
    
    figure;
    plot(xi,pdfi);
    step = xi(2) - xi(1);
    
    std = sigma;
    idxf(1) = min( find(xi>= 0.5*std) );
    idxi(1) = max( find(xi<=-0.5*std) );
    idxf(2) = min( find(xi>=   std) );
    idxi(2) = max( find(xi<=  -std) );
    idxf(3) = min( find(xi>= 2*std) );
    idxi(3) = max( find(xi<=-2*std) );
    idxf(4) = min( find(xi>= 3*std) );
    idxi(4) = max( find(xi<=-3*std) );
    
    for i= 1:4
        cums(i) = sum(pdfi(idxi(i):idxf(i))*step);
    end
    disp(cums)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
