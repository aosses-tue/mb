function [onset, dur, t_sil_aft, t_total_duration] = Create_noise_dau1996_default(nTag)
% function [onset, dur, t_sil_aft, t_total_duration] = Create_noise_dau1996_default(nTag)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       nTag = 3; % 600-ms white noise, no onset; no offset
%       [onset, dur, t_sil_aft, t_total_duration] = Create_noise_dau1996_default(nTag);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 23/10/2014 % Update this date manually
% Last use on   : 23/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
switch nTag
    case 1
        
        onset       = 100e-3;
        dur         = 300e-3;
        t_sil_aft   = 200e-3;
        t_total_duration = onset + dur + t_sil_aft; 
        
    case 2
        
        onset       = 100e-3;
        dur         = 300e-3;
        t_sil_aft   = 200e-3;
        t_total_duration = onset + dur + t_sil_aft; 
        
    case 3
        
        onset       = 0e-3;
        dur         = 600e-3;
        t_sil_aft   = 0e-3;
        t_total_duration = onset + dur + t_sil_aft; 
        
    case 10
        
        onset       = 0e-3;
        dur         = 200e-3;
        t_sil_aft   = 400e-3;
        t_total_duration = onset + dur + t_sil_aft; 
            
    case 20 % same as 10
        
        onset       = 0e-3;
        dur         = 200e-3;
        t_sil_aft   = 400e-3;
        t_total_duration = onset + dur + t_sil_aft; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
