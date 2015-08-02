function wavfilename = Get_LISTf_wav_filename(sentence,info)
% function wavfilename = Get_LISTf_wav_filename(sentence,info)
%
% Programmed by Alejandro Osses,  ExpORL, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavfilename = [];
pa          = getpaths;
tmp_char    = strsplit(sentence,' ');

listpath    = [pa.local_speechmaterials info.uri]; % x-Drive, ExpORL
txt_file    = 'LIST_zinnen.txt';
fid         = fopen([listpath txt_file]);

sentnr      = 0;
bMatched    = 0;

while bMatched == 0
    tline = fgetl(fid);
        
    tokens = regexp(tline, '([0-9]+). (.*\D)', 'tokens');
    if (strcmp(tokens{1}{2},sentence) )
        % sentnr = sentnr + 1;
        bMatched = 1;
    end
    sentnr=sentnr+1;
end

wavfilename = ['wdz' num2str(sentnr)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end