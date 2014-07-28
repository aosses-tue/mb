function note = num2note(num)
% function note = num2note(num)
%
% AOV added 0 for 'B'

switch num
    case 0
        note = 'B';
    case 1
        note = 'C';
    case 2
        note = 'C#';
    case 3
        note = 'D';
    case 4 
        note = 'D#';
    case 5 
        note = 'E';
    case 6
        note = 'F';
    case 7
        note = 'F#';
    case 8
        note = 'G';
    case 9 
        note = 'G#';
    case 10
        note = 'A';
    case 11
        note = 'A#';
    case 12
        note = 'B';
    otherwise
        error('Inproper number choosen');
end