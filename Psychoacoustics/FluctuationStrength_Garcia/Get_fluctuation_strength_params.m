function params = Get_fluctuation_strength_params(N,fs,dataset)
% function params = Get_fluctuation_strength_params(N,fs,dataset)
% 
% Creates the file 'params.mat' (if the file already exists it is deleted
% first) containing all the required parameters for the fluctuation
% strength model.
% 
%   - gzi was deleted on 07/01/2015
%   - Ndataset = 0; is the final fitting, as presented in the thesis
% 
% Author: Rodrigo Garcia
% Original name: Create_params (renamed when copied from Rodrigo's repository)
% Modified by: Alejandro Osses V.
% Last modified: 07/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dataset = 0; % 0 = Approved version
end
if nargin < 2
    fs = 44100;
end
if nargin < 1
    N = 12*fs; % 12 seconds
end

params         = struct;
params.fs      = fs;
params.N       = N;
params.Chno    = 47;
params.debug   = 'none';
switch dataset
    case {0,99}
        params.window_type = 'cosine';
    case 1
        params.window_type = 'blackman';
end
params.dataset = dataset;

params.Hweight = Get_Hweight_fluctuation(fs);
params.gzi     = Calculate_gzi(params.Chno,dataset);

switch dataset
    case 0
        % params.p_g     = 2;
        % params.p_m     = 0.25;
        % params.p_k     = 0.375;
        % params.cal = 0.1135; % New calibration factor by AO on 5/02/2016
        params.p_g     = 1;    % gzi = 1
        params.p_m     = 1.25;  % 1.25    1.5     2       2       2
        params.p_k     = 0.35;  % 0.4     0.4     0.4     0.5     2
        params.cal     = 0.2016;  % 0.2126  0.2477  0.335   0.3377  0.3621 New calibration factor by AO on 5/02/2016
    case 1
        params.p_g     = 2;
        params.p_m     = 0.25;
        params.p_k     = 0.375;
        params.cal     = 0.15;
    otherwise
        disp('Calibrating...assign your own values')
end
