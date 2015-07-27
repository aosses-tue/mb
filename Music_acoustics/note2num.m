function num = note2num(note)
% function num = note2num(note)
%
% Update for recognising Gsh and G#
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch note 
    case 'C'
        num = 1;
    case {'Csh','C#'}
        num = 2;
    case 'Db'
        num = 2;
    case 'D'
        num = 3;
    case {'Dsh','D#'}
        num = 4;
    case 'Eb'
        num = 4;
    case 'E'
        num = 5;
    case 'F'
        num = 6;
    case {'Fsh','F#'}
        num = 7;
    case 'Gb'
        num = 7;
    case 'G'
        num = 8;
    case {'Gsh','G#'}
        num = 9;
    case 'Ab'
        num = 9;    
    case 'A'
        num = 10;
    case {'Ash','A#'}
        num = 11;
    case 'Bb'
        num = 11;
    case 'B'
        num = 12;
    otherwise
        error('Inproper note expression');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%