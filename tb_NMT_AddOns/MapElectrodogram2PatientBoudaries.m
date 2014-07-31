function [v, txt_info, info] = MapElectrodogram2PatientBoudaries(v,p, FullScale)
% function [v, txt_info, info] = MapElectrodogram2PatientBoudaries(v,p, FullScale)
%
% p:    Structure containing information about Patient (T and C Levels)
% v:    Structure containing electrodogram with Current levels ranging from 
%       T and C Levels (constrained between 0 and 255 Clinical Units)
% txt_info: Output (text format) showing the parameters used in the scaling
%       (FAT, T- and C-Levels of the lowerst and highest frequency band
%
% Script validated to be used scaling electrodogramms using the Cochlear
% Ltd. nomenclature for FATs of 22 and 21 electrodes.
% 
% The resulting v-structure contains a new mapping per electrode in order 
% to get a current unit ranging from 0 and 255*FullScale where 0 means a CU 
% equal to T-Level and 255*FullScale  means a CU equal to C-Level. This new 
% structure is useful just for ploting ends (not for energy summation)
%
% Programmed by Alejandro Osses, ExpORL, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    FullScale = 0.8; % 80[%]
end

info.i_max = 22;
info.i_min = info.i_max+ 1 - p.CIC.NumChannels;

NumElectrode = info.i_max:-1:info.i_min;

for i = 1:length(NumElectrode)
    
    count = find(v.electrodes == NumElectrode(i));
    
    MinTL(i) = p.STIM.TLevels(i);
    MaxCL(i) = p.STIM.CLevels(i);
    MaxDR(i) = MaxCL(i)-MinTL(i);
    v.current_levels(count) = (v.current_levels(count) - MinTL(i))/MaxDR(i) * (256*FullScale);
    
end

txt_info = {['Electrodes: ' num2str(info.i_max) '  -   ' num2str(info.i_min) ' (FAT = ' num2str(length(NumElectrode)) ' channels)']; ...
            ['CLevels   : ' num2str(p.STIM.CLevels(1)) ' - ' num2str(p.STIM.CLevels(length(NumElectrode)))]; ...
            ['TLevels   : ' num2str(p.STIM.TLevels(1)) ' - ' num2str(p.STIM.TLevels(length(NumElectrode)))]};
info.NumElectrode = NumElectrode;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([mfilename '.m: Electrodogram scaled'])
disp(txt_info)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end