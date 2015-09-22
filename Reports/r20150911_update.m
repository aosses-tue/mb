function r20150911_update
% function r20150911_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 08/09/2015
% Last update on: 08/09/2015 
% Last use on   : 08/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bDoCrossCorrelation = 1;

if bDoCrossCorrelation
    N = 10000;
    M = 1;
    x = ones(N,M);
    
    fs = 44100;
    [xx xx xnorm] = Normalise_signal(x,fs); % fs not used for xnorm
    
    sigma1 = 0.8;
    x1 = Add_gaussian_noise(x,0,sigma1);
    
    sigma2 = 1.6;
    x2 = Add_gaussian_noise(x,0,sigma2);
    
    sigma3 = 2.4;
    x3 = Add_gaussian_noise(x,0,sigma3);
    
    sigma4 = 50;
    x4 = Add_gaussian_noise(x,0,sigma4);
    X = [x x1 x2 x3 x4];
    
    std(X)
    
    Method = 1;
    if Method == 1
        cc(1) = optimaldetector(x,xnorm);
        cc(2) = optimaldetector(x1,xnorm);
        cc(3) = optimaldetector(x2,xnorm);
        cc(4) = optimaldetector(x3,xnorm);
        cc(5) = optimaldetector(x4,xnorm);
        
        xnormn = Add_gaussian_noise(xnorm,0,0.8);
        [ccn(1) xcorre] = optimaldetector(x,xnormn);
        
        bGaussian(1,1) = Is_normal_distributed(x);
        bGaussian(1,2) = Is_normal_distributed(xnormn);
        bGaussian(1,3) = Is_normal_distributed(xcorre);
        
        xnormn = Add_gaussian_noise(xnorm,0,0.8);
        [ccn(2) xcorre] = optimaldetector(x1,xnormn);
        
        bGaussian(2,1) = Is_normal_distributed(x1);
        bGaussian(2,2) = Is_normal_distributed(xnormn);
        bGaussian(2,3) = Is_normal_distributed(xcorre);
        
        xnormn = Add_gaussian_noise(xnorm,0,0.8);
        ccn(3) = optimaldetector(x2,xnormn);
        xnormn = Add_gaussian_noise(xnorm,0,0.8);
        ccn(4) = optimaldetector(x3,xnormn);
        xnormn = Add_gaussian_noise(xnorm,0,0.8);
        ccn(5) = optimaldetector(x4,xnormn);
    elseif Method == 2
        Method = 'coefficient-non-normalised';
        cc(1) = Correlation(x,xnorm,Method);
        cc(2) = Correlation(x1,xnorm,Method);
        cc(3) = Correlation(x2,xnorm,Method);
        cc(4) = Correlation(x3,xnorm,Method);
        cc(5) = Correlation(x4,xnorm,Method);
    elseif Method == 3
        Method = 'coefficient';
        cc(1) = Correlation(x,xnorm,Method);
        cc(2) = Correlation(x1,xnorm,Method);
        cc(3) = Correlation(x2,xnorm,Method);
        cc(4) = Correlation(x3,xnorm,Method);
        cc(5) = Correlation(x4,xnorm,Method);
    end

end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
