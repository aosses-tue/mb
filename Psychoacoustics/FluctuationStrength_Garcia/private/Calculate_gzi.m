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
        % g0 = [
        % 0       1
        % 125     1
        % 250     1
        % 500     1
        % 1000    1
        % 1500    1
        % 2000    1
        % 3000    1
        % 4000    1
        % 6000    1
        % 8000    1
        % 16000   1
        % ];

%     g0 = [  0  ,1  ,2.5,4.9,6.5,8,9,10  ,11,11.5,13,15  ,17.5 ,24;
%             0.5,0.9,1  ,1  ,1  ,1,1,1   , 1, 1  , 1, 0.9, 0.7 ,0.5];
    g0 = [  0,1,2.5,4.9,6.5,8,9,10  ,11,11.5,13,15  ,17.5 ,24;
            1,1,1  ,1  ,1  ,1,1,1   , 1, 1  , 1, 0.9, 0.7 ,0.5];
	g0 = transpose(g0);
    
    % for k = 1:47
    %     gzi(k)  = sqrt(interp1(gr(1,:)',gr(2,:)',k/2));
    % end
    gzi = interp1(g0(:,1),g0(:,2),(1:Chno)*Chstep);
    gzi(isnan(gzi)) = g0(end,2); % 0

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
    gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
    gzi(isnan(gzi)) = g0(end,2); % 0
    error('Check')
    case 90
        % Looking at values AM fc above 1000 Hz 
        %           1-([1.26 1.37 1.39 0.86]-1.2)/1.2
     g0 = [
        0       1
        1000    1
        8000    0.8
        16000   0.7
    ];
    % gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
    % gzi(isnan(gzi)) = g0(end,2); % 0
    error('Check')
end   
    
end
