function     out = nsb_ReadFAT(FATn)
% Routine to read the FAT table from xml file. Function returns the channel
% edge frequencies (in Hz) in a vector. FATn specifies the FAT table number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Company Confidential. Copyright(c) Cochlear Limited 2004-2005.
%  Nucleus Matlab Blockset (NMB).  For use with NATB S/W and H/W
%
%  $Archive: $
%  $Revision: #1 $
%  $Date: 2011/05/02 $
%  Authors:     John Heasman, Adam Hersbach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % read FAT tables from the FAT.xml
    xmlstr  = fileread('FAT.xml'); 
    V       = xml_parseany(xmlstr);
    
    table_sel   = num2str(FATn);
    if (str2double(table_sel) < 1) || (str2double(table_sel) >22)           % yanyou
        if str2double(table_sel) > 0
            errordlg('FATs has to be between 2 to 22','FAT error','modal')  % yanyou
            return;
        end
    end
    
    % extract FAT frequency band information
    for i = 1:length(V.FATTable)
        if (strcmp(V.FATTable{i}.name{1}.CONTENT, table_sel)) == 1
            break;
        end
    end  
    
    freq_band = V.FATTable{str2double(table_sel)}.FATBand;
    
    for i = 1:length(freq_band)
        char_freqs(i) = (str2double(freq_band{i}.LowerFrequency{1}.CONTENT) + ...
                                    str2double(freq_band{i}.UpperFrequency{1}.CONTENT))/2; 
        crossover_freqs(i) = str2double(freq_band{i}.LowerFrequency{1}.CONTENT)';
    end    
    crossover_freqs(i+1) = str2double(freq_band{i}.UpperFrequency{1}.CONTENT)';
    
    chan_edges  = [crossover_freqs(1:end-1)', crossover_freqs(2:end)'];
    
    out         = crossover_freqs;
    
    
