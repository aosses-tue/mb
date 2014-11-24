function q = Decode_embedded_protocol(f)

% Decode_embedded_protocol: Decode an embedded protocol sequence.
%
% q = Decode_embedded_protocol(f)
%
% f:	Embedded protocol frame sequence (fields listed below).
%
% q:    Pulse sequence (fields listed below).
%
%             
% Protocol frame sequence fields    Pulse sequence fields                    
% ------------------------------    ---------------------      
% es                                electrodes                        
% ms                                modes                             
% as                                current_levels                    
% ws (RF cycles)                    phase_widths (microseconds)       
% gs (RF cycles)                    phase_gaps   (microseconds)       
% ts (RF cycles)                    periods      (microseconds)       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROTOCOL = Embedded_protocol_params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode & Mode:

[q.electrodes, q.modes] = Decode_embedded_electrodes(f.es, f.ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current level

q.current_levels = f.as;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase width

phase_width_cycles = f.ws + PROTOCOL.PHASE_WIDTH_BASE_CYCLES;
q.phase_widths = phase_width_cycles / 5;		% 5 MHz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase gap

phase_gap_cycles = f.gs + PROTOCOL.PHASE_GAP_BASE_CYCLES;
q.phase_gaps = phase_gap_cycles / 5;			% 5 MHz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Period

q.periods = f.ts / 5;							% 5 MHz.
