function Save_sequence(seq, base_name)

% Save_sequence: Saves a sequence to a file.
%
% function Save_sequence(seq, base_name)
%
% Inputs:
% seq:          Sequence (channel-magnitude or frame sequence).
% base_name:    Base file name.
%
% The sequence is saved to disk dependent upon the sequence type, which is
% based upon the structure of the sequence.  The supported formats and their
% file extensions are as follows:
%   .qic    Channel/magnitude sequence.
%   .quf    SP5 format sequence.
% A sequence file is a little-endian binary file of 16-bit unsigned integers,
% consisting of the following command tokens:
% Frame token:
%   E, R, CL, PWCNT, PWSEL
% Channel/magnitude token:
%   6 MSB:      channel, range (1 - 23)
%   10 LSBs:    magnitude
% Period token:
%   1 MSB:      1
%   15 LSBs:    period in 5 MHz cycles
%
% At present only CI24M 5 MHz embedded protocol uniform rate frame sequences
% are supported (.quf).
%
% Examples:
%   Save_sequence(seq, 'foo');      % creates foo.quf
%   Save_sequence(seq);             % creates seq.quf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson, Colin Irwin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SEQ_INFO = Sequence_tokens;

file_type = [];

% Determine the sequence type:

if (isfield(seq, 'channels') & ...
    isfield(seq, 'magnitudes'))
    
    % The sequence is a channel/magnitude (floating point) sequence.
    file_type = 'qic';
    image_id = SEQ_INFO.QIC_TYPE;
    
elseif (isfield(seq, 'electrodes')      & ...
        isfield(seq, 'modes')           & ...
        isfield(seq, 'current_levels')  & ...
        isfield(seq, 'phase_widths')    & ...
        isfield(seq, 'phase_gaps')      & ...
        isfield(seq, 'periods'))
        
    % The sequence is a frame sequence.
    file_type = 'quf';
    image_id = SEQ_INFO.QUF_TYPE;
    
else    
    error('Unsupported sequence type');
end

file_name = strcat(base_name, '.', file_type);

% Open the file in little-endian format ('l' flag):
fid = fopen(file_name,'w','l');
if (fid == -1)
    error('Cannot open file');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the sequence to file, appropriately for the type.
if (file_type == 'qic')
    for i = 1:(length(seq.channels)),
        fwrite(fid, seq.channels(i), 'uint16');
        fwrite(fid, seq.magnitudes(i), 'float32');
    end    
else

    [rfseq, idle] = Encode_embedded_protocol(seq);

    % Write sequence image type token:
    fwrite(fid, SEQ_INFO.TOKEN_IMAGE_TYPE, 'uint16');
    fwrite(fid, image_id, 'uint16');

    fwrite(fid, SEQ_INFO.TOKEN_PROTOCOL_CONFIG, 'uint16');
    fwrite(fid, 0, 'uint16');

    % Constant phase_gap:
    if length(seq.phase_gaps) > 1
        error('Only sequences with constant phase gap are supported');
    end
    phase_gap_token = Calc_phase_gap_token(rfseq, SEQ_INFO);
    fwrite(fid, phase_gap_token, 'uint32');

    % Constant period:
    if length(seq.periods) > 1
        error('Only uniform-rate sequences are supported');
    end
    period_token = Calc_period_token(rfseq, SEQ_INFO);
    fwrite(fid, period_token, 'uint16');

    % Idle frame:
    idle_frame_token = Calc_frame_token(idle, SEQ_INFO);
    fwrite(fid, idle_frame_token, 'uint32');

    % Pulses:
    frame_tokens = Calc_frame_token(rfseq, SEQ_INFO);
    fwrite(fid, frame_tokens, 'uint32');
    
    % End-of-sequence:
    % We really only need one 16-bit zero word,
    % but we write two words to make Read_sequence easier.
    fwrite(fid, 0, 'uint32');
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function period_token = Calc_period_token(fseq, SEQ_INFO)

if any(fseq.ts >= (2 ^ 14))       % count in bits 0-13, clock divider = 1 = '00'
    error('period is too large');
end

period_token  = SEQ_INFO.TOKEN_PERIOD + fseq.ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phase_gap_token = Calc_phase_gap_token(fseq, SEQ_INFO)

if any(fseq.gs >= (2 ^ 8))        % count in bits 0-7, clock divider = 1 = '00'
    error('phase gap is too large');
end

% 32-bit token: little-endian file has the least significant word written first:
phase_gap_token = SEQ_INFO.TOKEN_PHASE_GAP + bitshift(fseq.gs, 16);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame_token = Calc_frame_token(fseq, SEQ_INFO)

% DEF frame word 0: E, M, assume W_sel = 0:

frame_0 = SEQ_INFO.TOKEN_FRAME + bitshift(fseq.ms, 5) + fseq.es;

% DEF frame word 1: W_count, A:

if any(fseq.ws >= (2 ^ 8))
    error('phase width is too large');
end
frame_1 = bitshift(fseq.ws, 8) + fseq.as;

% 32-bit token: little-endian file has the least significant word written first:
frame_token = frame_0 + bitshift(frame_1, 16);
