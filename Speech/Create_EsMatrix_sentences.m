function Create_EsMatrix_sentences
% function Create_EsMatrix_sentences
%
% 1. Description:
%       Creates the sentences of the Spanish Matrix as recorded by AO. It 
%       pastes the name-verb-numeral-subject and adjective, which where 
%       previously recorded.
% 
% 2. Stand-alone example:
%       Create_EsMatrix_sentences;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 14/06/2015
% Last update on: 14/06/2015 % Update this date manually
% Last use on   : 14/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

pathmatrix = [Get_TUe_paths('db_speechmaterials') 'spanish' delim 'Matrix' delim];

pathmatrixout = [Get_TUe_paths('db_speechmaterials') 'spanish' delim 'Matrix' delim 'New-files' delim];
Mkdir(pathmatrixout);

srcpath = [pathmatrix '01-Nuendo-session' delim 'Matrix-ES-edit' delim];
subpath = {'01-Nombre'; '02-Verbo'; '03-Numeral'; '04-Sustantivo'; '05-Adjetivo'};

m = es_speech_material('Matrix');

[X N M] = size(m);

for in = 1:N
    for im = 1:M
        st = m{1,in,im};
        
        [f1 fs] = Wavread([srcpath subpath{1} delim st.file(1) '.wav']);
        f2      = Wavread([srcpath subpath{2} delim st.file(2) '.wav']);
        f3      = Wavread([srcpath subpath{3} delim st.file(3) '.wav']);
        f4      = Wavread([srcpath subpath{4} delim st.file(4) '.wav']);
        f5      = Wavread([srcpath subpath{5} delim st.file(5) '.wav']);
        
        x = [f1; f2; f3; f4; f5];
        
        Wavwrite(x,fs,[pathmatrix st.file]);
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
