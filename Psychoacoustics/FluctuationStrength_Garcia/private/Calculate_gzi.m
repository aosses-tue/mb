function gzi = Calculate_gzi(Chno,dataset)
% function gzi = Calculate_gzi(Chno,dataset)
%
% Returns gzi parameters using the specified number of channels.
% 
% Inputs:
% Chno: The number of channels.
% 
% Outputs:
% gzi: The gzi vector.
% 
% Author: Rodrigo Garcia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dataset = 0;
end
Chstep = 0.5;
    
switch dataset
    case {0,99}
        % Looking at values AM fc above 1000 Hz 
        %           1-([1.26 1.37 1.39 0.86]-1.2)/1.2
%         g0 = [
%         0       0.8
%         125     0.8
%         250     0.85
%         500     0.9
%         1000    0.9500 % intuitively up to here
%         1500    0.9
%         2000    0.8583    
%         3000    0.85
%         4000    0.8417    
%         6000    0.84 % 1
%         8000    0.84 % 1.2833
%         16000   0.84
%     ];
%         g0 = [
%         0       0.7727
%         125     0.7727
%         250     0.8463
%         500     0.9 %0.9220
%         1000    1     
%         1500    0.9220
%         2000    0.8587
%         3000    0.8463
%         4000    0.8339
%         6000    0.8314
%         8000    0.8314
%         16000   0.8314
%     ];
g0 = [
        0       1
        125     1
        250     1
        500     1 
        1000    1     
        1500    0.9220
        2000    0.8587
        3000    0.8463
        4000    0.8339
        6000    0.8314
        8000    0.8314
        16000   0.8314
        ]; 
    
    case 1
    g0 = [
        0       0
        125     1
        250     1
        500     1
        1000    1
        1500    1
        2000    1
        3000    1
        4000    1
        6000    1
        8000    1
        16000   0
    ];
end   
gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
gzi(isnan(gzi)) = 0;

% gzi = ones(1,Chno);
    
end
