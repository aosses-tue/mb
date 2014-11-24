function SEQ_INFO = Sequence_tokens
% Sequence_tokens: Returns a struct containing sequence token constants
% function SEQ_INFO = Sequence_tokens

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image types:

SEQ_INFO.PERIOD_TYPE            =  256;
SEQ_INFO.FRAME_TYPE             =  512;
SEQ_INFO.INT_CHAN_MAG_TYPE      = 1024;
SEQ_INFO.FLOAT_CHAN_MAG_TYPE    = 2048;

SEQ_INFO.SEQ_IMAGE_VERSION      =    2;

% Uniform rate, frame tokens.
SEQ_INFO.QUF_TYPE               = SEQ_INFO.FRAME_TYPE + SEQ_INFO.SEQ_IMAGE_VERSION;

% Uniform rate, channel / magnitude (integer - SPrint compatible) tokens.
SEQ_INFO.QUC_TYPE               = SEQ_INFO.INT_CHAN_MAG_TYPE + SEQ_INFO.SEQ_IMAGE_VERSION;

% Uniform rate, channel / magnitude (floating point) tokens.
SEQ_INFO.QIC_TYPE               = SEQ_INFO.FLOAT_CHAN_MAG_TYPE + SEQ_INFO.SEQ_IMAGE_VERSION;

% Token IDs:

SEQ_INFO.TOKEN_IMAGE_TYPE           =  1;
SEQ_INFO.TOKEN_PROTOCOL_CONFIG      =  2;
SEQ_INFO.TOKEN_PHASE_GAP            = 16;

SEQ_INFO.TELEMETRY_NBITS            =  8;
SEQ_INFO.MAGNITUDE_NBITS            = 10;
SEQ_INFO.FRAME_NBITS                = 13;
SEQ_INFO.PERIOD_NBITS               = 15;

SEQ_INFO.TOKEN_TELEMETRY            = bitshift(1, SEQ_INFO.TELEMETRY_NBITS);
SEQ_INFO.TOKEN_CHANNEL_MAGNITUDE    = bitshift(1, SEQ_INFO.MAGNITUDE_NBITS);
SEQ_INFO.TOKEN_FRAME                = bitshift(3, SEQ_INFO.FRAME_NBITS);
SEQ_INFO.TOKEN_PERIOD               = bitshift(1, SEQ_INFO.PERIOD_NBITS);


