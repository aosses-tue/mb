function [note, GluedNote] = freq2note(freq)
% function [note, GluedNote] = freq2note(freq)
%
% Programmed by Matthias, updated by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return;
end

times_divided = 0;

if freq > 440
    % Formula does not work for frequencies above 440 Hz, then we divide by
    % 2 until the frequency lies below 440 Hz. The octave number is then
    % corrected summing 1 octave per each time we divided by 2
    while freq > 440
        freq = freq/2;
        times_divided = times_divided + 1;
    end
end

note.octave     = 4;
note.note       = 'A'; % 440 Hz, A4
Num_octaves     = log10(freq / 440) / log10(2);
note.octave     = note.octave - floor(abs(Num_octaves)) + times_divided;

SemitonesBelowA = 12-round(mod( Num_octaves,1 )*12);

switch SemitonesBelowA
    case 0
        note.note = 'A';
    case 1
        note.note = 'Gsh';
    case 2
        note.note = 'G';
    case 3
        note.note = 'Fsh';
    case 4
        note.note = 'F';
    case 5
        note.note = 'E';
    case 6
        note.note = 'Dsh';
    case 7
        note.note = 'D';
    case 8
        note.note = 'Csh';
    case 9
        note.note = 'C';
    case 10
        note.note = 'B';
    case 11
        note.note = 'Ash';
end

GluedNote = [note.note num2str(note.octave)];
