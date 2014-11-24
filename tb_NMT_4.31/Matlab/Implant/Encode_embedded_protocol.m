function [f, idle] = Encode_embedded_protocol(q)

% Encode_embedded_protocol: Encode a sequence into embedded protocol.
%
% [f, idle] = Encode_embedded_protocol(q)
%
% q:    Pulse sequence (fields listed below).
%
% f:    Embedded protocol frame sequence (fields listed below).
% idle: Idle frame inferred from the pulse sequence.
%       In MP1+2, a special idle frame can be used: (E,M,A) = (24,25,0).
%       In MP1 or MP2, there is no special idle mode, instead
%       the most apical electrode is used with CL = 0.
%
%
% Pulse sequence fields             Protocol frame sequence fields
% ---------------------------       ------------------------------
% electrodes                        es
% modes                             ms
% current_levels                    as
% phase_widths (microseconds)       ws (RF cycles)
% phase_gaps   (microseconds)       gs (RF cycles)
% periods      (microseconds)       ts (RF cycles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.implant = Ensure_CIC3_params;
Check_sequence(p, q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrode & Mode:

special_idles = find(q.electrodes == 0);

% Check the modes:
if length(q.modes) > 1
    error('Only constant mode is supported');
end

MODE = Implant_modes;

switch (q.modes)

    case MODE.MP1
        f.es = q.electrodes;
        f.ms = 24;
        idle.es = max(f.es);    % most apical electrode used.
        idle.ms = f.ms;
        if any(special_idles)
            error('Special idle frames not available in MP1');
        end
        
    case MODE.MP2
        f.es = q.electrodes;
        f.ms = 25;
        idle.es = max(f.es);    % most apical electrode used.
        idle.ms = f.ms;
        if any(special_idles)
            error('Special idle frames not available in MP2');
        end
        
    case MODE.MP1+2
        f.es = q.electrodes;
        f.ms = 30;
        idle.es = 24;
        idle.ms = 25;
        if any(special_idles)
			num_pulses = length(q.electrodes);
            f.ms = repmat(30, num_pulses, 1);
            f.es(special_idles) = idle.es;
            f.ms(special_idles) = idle.ms;
        end
        
    otherwise
        error('Unsupported mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current level

f.as = q.current_levels;
idle.as = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase width

phase_width_cycles = round(q.phase_widths * 5);     % 5 MHz
f.ws = phase_width_cycles - p.implant.PHASE_WIDTH_BASE_CYCLES;
idle.ws = min(f.ws);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase gap

if ~all(q.phase_gaps >= p.implant.MIN_PHASE_GAP_us)
    error('Phase gap is too small');
end
if ~all(q.phase_gaps <= p.implant.MAX_PHASE_GAP_us)
    error('Phase gap is too large');
end

phase_gap_cycles = round(q.phase_gaps * 5);         % 5 MHz.
f.gs = phase_gap_cycles - p.implant.PHASE_GAP_BASE_CYCLES;
idle.gs = min(f.gs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Period

if ~all(q.periods >= p.implant.MIN_PERIOD_us)
    error('Period is too small');
end
if ~all(q.periods <= p.implant.MAX_PERIOD_us)
    error('Period is too large');
end

f.ts = round(q.periods * 5);                        % 5 MHz.
idle.ts = min(f.ts);
