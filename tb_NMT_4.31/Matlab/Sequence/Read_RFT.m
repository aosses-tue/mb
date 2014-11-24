function seq = Read_RFT(filename)

% Read_RFT: Read an RFT file into a pulse sequence struct.
% An RFT file is a text file created by RxFrames.exe,
% containing frame information captured by the PCI.
% Note that phase 2 width is ignored.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename,'r');
if (fid == -1)
    error('Cannot open file');
end

pulses              = fscanf(fid, '%d %d %d %f %f %f %f', [7,inf])';
es                  = pulses(:,1);
ms                  = pulses(:,2);
[seq.electrodes, seq.modes] = Decode_embedded_electrodes(es, ms);
seq.current_levels  = pulses(:,3);
seq.phase_widths    = pulses(:,4);
seq.phase_gaps      = pulses(:,5);
seq.periods         = pulses(:,7);
