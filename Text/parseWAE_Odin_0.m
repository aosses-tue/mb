function M = parseWAE_Odin_0
% function M = parseWAE_Odin_0
%
% 1. Description:
%
% 2. Stand-alone example:
%     M = parseWAE_Odin_0;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'D:\Documenten-TUe\06-Lectures-TUe\BEP-2015-2016\03-Odin\Experiment\webaudioevaluationtool-45df85a336d4\saves\';

files = {'save-AO.xml','save-AO-proc.xml','test-1.xml'}
% pairs = [0 1; 0 2; 0 3; 0 4; 2 1; 3 1; 4 1; 2 3; 2 4; 3 4]; % put the pairs here, manually

for i = 1:length(files)
    resultfile = [dir files{i}];
    
    M = parseWAE_Odin_1(resultfile);
    
    if i == 1
        Mtot = M;
    else
        Mtot = Mtot + M;
    end
    
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end