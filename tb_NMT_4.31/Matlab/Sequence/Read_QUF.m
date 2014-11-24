function seq = Read_QUF(file_name)

% Read_QUF: Read a QUF file into a pulse sequence struct.
% A QUF file is a little-endian binary file of uint16,
% consisting of SP5 NIC tokens.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SEQ_INFO = Sequence_tokens;
SEQ_INFO.PROTOCOL = Embedded_protocol_params;
image_type = SEQ_INFO.QUF_TYPE;

fid = fopen(file_name,'r','l'); % Little-endian (Intel format)
if (fid == -1)
    error('Cannot open file');
end

% Create the fields in the order we like to display them:

seq.electrodes      = [];
seq.modes           = [];
seq.current_levels  = [];
seq.phase_widths    = [];
seq.phase_gaps      = [];
seq.periods         = [];

% Read sequence image type token from binary file:
type_token  = fread(fid, 1, 'uint16');
type_value  = fread(fid, 1, 'uint16');
if ((type_token ~= SEQ_INFO.TOKEN_IMAGE_TYPE) | (type_value ~= image_type))
    type_token
    type_value
    error('File has incorrect type_token');
end

protocol_token = fread(fid, 1, 'uint16');
protocol_value = fread(fid, 1, 'uint16');
if ((protocol_token ~= SEQ_INFO.TOKEN_PROTOCOL_CONFIG) | (protocol_value ~= 0))
    protocol_token
    protocol_value
    error('File has incorrect protocol_token');
end

% Constant phase_gap:
phase_gap_token = fread(fid, 1, 'uint32');
fseq.gs = Decode_phase_gap_token(phase_gap_token, SEQ_INFO);

% Constant period:
period_token = fread(fid, 1, 'uint16');
fseq.ts = Decode_period_token(period_token, SEQ_INFO);

% Idle frame:
% Just check that the CL is zero.
idle_frame_token = fread(fid, 1, 'uint32');
[idle_es, idle_ms, idle_as, idle_ws] = Decode_SP5_frame_token(idle_frame_token);
if (idle_as ~= 0)
    idle_as
    error('Idle CL > 0')
end

% Pulses:
frame_tokens = fread(fid, inf, 'uint32');
% The last 32-bit word should be all zero (denoting End-of-sequence):
if (frame_tokens(end) ~= 0)
    error('Missing End-of-sequence token');
end
frame_tokens(end) = []; % discard End-of-sequence token.
[fseq.es, fseq.ms, fseq.as, fseq.ws, id] = Decode_SP5_frame_token(frame_tokens);

if ~all(id == 3)
    error('Bad frame token');
end

fclose(fid);

seq = Decode_embedded_protocol(fseq);

if all(seq.phase_widths == seq.phase_widths(1))
    seq.phase_widths = seq.phase_widths(1);
end
if all(seq.modes == seq.modes(1))
    seq.modes = seq.modes(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-functions used in Read_QUF:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function period_cycles = Decode_period_token(token, SEQ_INFO)

if (bitand(token, SEQ_INFO.TOKEN_PERIOD) ~= SEQ_INFO.TOKEN_PERIOD)
    error('Bad period token');
end

period_cycles = token - SEQ_INFO.TOKEN_PERIOD;
if any(period_cycles >= (2 ^ 14))       % count in bits 0-13, clock divider = 1 = '00'
    error('period is too large');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function protocol_phase_gap_cycles = Decode_phase_gap_token(token, SEQ_INFO)

if (bitand(token, SEQ_INFO.TOKEN_PHASE_GAP) ~= SEQ_INFO.TOKEN_PHASE_GAP)
    error('Bad phase gap token');
end

protocol_phase_gap_cycles = bitshift(token, -16);   
if any(protocol_phase_gap_cycles >= (2 ^ 8))        % count in bits 0-7, clock divider = 1 = '00'
    error('phase gap is too large');
end
