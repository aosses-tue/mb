function wavfilename = Get_VlMatrix_wav_filename(sentence)
% function wavfilename = Get_VlMatrix_wav_filename(sentence)
%
% Example:
%       sentence = 'David draagt drie blauwe bedden';
%       wavfilename = Get_VlMatrix_wav_filename(sentence);
%       % The answer is going to be: '00110'
% 
% Programmed by Alejandro Osses,  ExpORL, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavfilename = [];

tmp_char = strsplit(sentence,' ');

Matrix{1} = {'David'; 'Ellen'; 'Emma'; 'Jacob';'Jeroen';'Johan';'Lucas';'Sara';'Sofie';'Thomas'};
Matrix{2} = {'draagt';'heeft';'kiest';'koopt';'krijgt';'leent';'maakt';'wint';'ziet';'zoekt'};
Matrix{3} = {'twee';'drie';'vier';'vijf';'zes';'acht';'tien';'elf';'twaalf';'veel'};
Matrix{4} = {'beige';'blauwe';'bruine';'gele';'grijze';'groene';'paarse';'rode';'witte';'zwarte'};
Matrix{5} = {'bedden';'boten';'doeken';'dozen';'fietsen';'jassen';'kousen';'manden';'pennen';'ringen'};

for i = 1:length(tmp_char)
    for j=1:10
        bFound = strcmp(Matrix{i}(j,1),tmp_char{i});
        if bFound == 1
            wavfilename = [wavfilename num2str(j-1)];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
