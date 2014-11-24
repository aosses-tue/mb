function [es, ms, as, ws, id] = Decode_SP5_frame_token(token)

% Decode_SP5_frame_token: Decode an SP5 frame token.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 32-bit token: little-endian format has the least significant word written first:

es      = Extract_bit_field(token,  0, 5);
ms      = Extract_bit_field(token,  5, 5);
psel    = Extract_bit_field(token, 10, 3);
id      = Extract_bit_field(token, 13, 3);
as      = Extract_bit_field(token, 16, 8);
ws      = Extract_bit_field(token, 24, 8);

if any(psel ~= 0)
    error('Only phase width multiplier 0 supported');
end
