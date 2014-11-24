function p = Embedded_protocol_params(p)
% Embedded_protocol_params: Calculates parameter struct for embedded protocol
% function p = Embedded_protocol_params(p)
% The struct contains various embedded protocol constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The RF waveform received at the implant is a few cycles longer
% than the transmitted RF waveform, due to ringing.
% This effect extends the current pulse.
% (It is not included in RF frame width calculations)

p.RINGING_CYCLES = 3;

% Number of cells per token.

p.TOKEN_CELLS = 6;

% Number of tokens per phase.

p.PHASE_TOKENS = 5;

% Number of cells per phase (before phase extender).
% The embedded protocol requires at least one cell to follow the parity token.
% The SP5 DEF inserts this extra cell automatically, so here
% the Phase Extender is defined to start after this extra cell.

p.PHASE_CELLS = p.PHASE_TOKENS * p.TOKEN_CELLS + 1;

% Number of cells (& cycles) before current pulse commences.
% The stimulus current does not start until one cell into the second token:

p.PHASE_DELAY_CELLS  = p.TOKEN_CELLS + 1;
p.PHASE_DELAY_CYCLES = p.PHASE_DELAY_CELLS * 5;

% Phase gap calculation offset:

p.PHASE_GAP_BASE_CYCLES = p.PHASE_DELAY_CYCLES - p.RINGING_CYCLES;

% Number of cells in current pulse before phase extender is added.
% The stimulus current phase width is the phase extender plus the initial width
% due to the tokens. The phase_delay_cells must be subtracted.

p.PHASE_INITIAL_CELLS  = p.PHASE_CELLS - p.PHASE_DELAY_CELLS;
p.PHASE_INITIAL_CYCLES = p.PHASE_INITIAL_CELLS * 5;

% Phase width calculation offset:

p.PHASE_WIDTH_BASE_CYCLES = (p.PHASE_INITIAL_CELLS * 5) + p.RINGING_CYCLES;

% Minimum period in microseconds:

p.MIN_PERIOD_us = 69.4;
p.MAX_PERIOD_us = 3276.6;

% Minimum phase gap in microseconds:

p.MIN_PHASE_GAP_us = 8;
p.MAX_PHASE_GAP_us = 57.4;
