function [Note, gluedNote, freqNote] = sumSemitones2note(note, Semitones)
% function [Note, gluedNote, freqNote] = sumSemitones2note(note, Semitones)
% 
% Inputs:
%   note        - structure with the fields note and octave
%   Semitones   - Semitones to be added to note. 
% 
% Example:
%   note.note   = 'E'; 
%   note.octave = '3'; 
%   Semitones   = 1;
%   [Note, gluedNote, freqNote] = sumSemitones2note(note,Semitones);
%
%   >> Expected result: Note.note = 'F';
%                       Note.octave = 2;    
%                       gluedNote = 'F3'; %(1 semitone above E3)
%                       freqNote = 175; % value in Hz
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reference.note  = note.note;
if ischar( note.octave ) == 1
    reference.octave = str2num(note.octave);
else
    reference.octave = note.octave;
end

freq_ref        = note2freq(reference, 1);

ref_midi        = note2num(reference.note) + reference.octave*12;

new_note_midi   = ref_midi + Semitones;
a               = 2^(1/12); % Amount of octaves between 2 consecutive semitones

ref_midi    = note2num(reference.note) + reference.octave*12;

halfsteps   = Semitones;
freqNote    = round( freq_ref*a^halfsteps );

[Note, gluedNote] = freq2note( freqNote );
