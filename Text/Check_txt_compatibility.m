function Check_txt_compatibility(filename)
% function Check_txt_compatibility(filename)
%
%
% % Example:
%       filename = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-generated-20140527-0817\prueba-4.txt';
%       Check_txt_compatibility(filename);
% 
% Programmed by Alejandro Osses
% Created on     : 04/06/2014
% Last updated on: 04/06/2014
% Last used on   : 04/06/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
%     [exp2 exp1] = uigetfile('*.txt','Select a txt file containing numeric data: ');
%     filename = [exp1 exp2];
    filename = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data\mode-1-v_2.txt';
end

literal         = '^'; % Non compatible characters
count           = 0;
nline           = 0;
disp(filename);

fid             = fopen(filename,'r');
text2print{1}   = '';
i_count         = 1;

tic
while ~feof(fid)
    nline       = nline+1;
    tline       = fgetl(fid);
    matches     = strfind(tline, literal);
    
    if length(tline)~=0
        if ~isempty(matches); % 'then it is compatible'

            aux = sprintf([text2print{i_count} '%s\n'],'0.00');
            count = count+1;

        else

            aux = sprintf([text2print{i_count} '%s\n'],tline);
            
        end
        
        text2print{i_count} = aux;
            
        if mod(nline,50000) == 0
            toc
            i_count = i_count+1;
            text2print{i_count} = [];
            tic
        end
        
    else
        break;
    end
    
    if mod(nline,10000)==0
        disp( num2str(nline) )
    end

end

fclose(fid);
toc
if count > 0
    % disp(text2print)
    Comp_filename = [Delete_extension(filename,'txt') '-c.txt'];
    fid = fopen(Comp_filename,'w'); % Creates and write data
    
    for i = 1:length(text2print)
        fprintf(fid,'%s',text2print{i});
    end
    
    fclose(fid);
    
    if fid > 0
        disp([Comp_filename ' created successfully...']);
    end
end

disp( [num2str(count) ' changes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end