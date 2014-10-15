function r20141009_KUL_IJA
% function r20141009_KUL_IJA
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 09/10/2014
% Last update on: 09/10/2014 % Update this date manually
% Last use on   : 09/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.bSave = 1;
options.bDoLISTwhite    = 1;
options.bDoLISTssn      = 1;
options.bAssess         = 0;
options.bDoPlotSNR      = 0;
options.bDoSpeech       = 0;
options.bDoPlotF0       = 0;
options.bDoAdditionalAnalysis = 0;

%%
options.bDoPDA = 1;
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure-IJA-F0ref-clean' delim];
options.bCleanSpeech    = 1;
r20140930_FIA_KUL(options);

%%
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure-IJA-F0ref-CP810' delim];
options.bCleanSpeech    = 0;
r20140930_FIA_KUL(options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
