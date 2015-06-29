function params = FluctuationStrength_Garcia_getparams(N)
% function params = FluctuationStrength_Garcia_getparams(N)
% 
% Returns a struct with all the required parameters for the fluctuation
% strength model.
% 
% Outputs:
% params: Parameters struct
%
% Author: Rodrigo Garcia
% Original name: Get_params (renamed when copied from Rodrigo's repository)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

PARAMS_FILE = sprintf('params-%.0f.mat',N);

if ~exist(PARAMS_FILE,'file')
    FluctuationStrength_Garcia_createparams(N); 
end

params = load(PARAMS_FILE);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
