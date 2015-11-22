function freq = note2freq(test,bRounded)
% function freq = note2freq(test, bRounded)
% 
% 1. Description:
%       test is a structure with two fields:
%           test.note (char)
%           test.octave (int)
%       freq is the frequency in Hz related to the tone defined by test. 
%       For example, if test.note is 'A', and test.octave is 4 then freq 
%       will be 440 (Hz)
%
% 2. Stand-alone example:
%   test.note = 'A'; 
%   test.octave = 4; 
%   note2freq(test);
% 
% Programmed by Matthias M.
% Original location: ../matthias/src/Matlab/Misc/note2freq.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <2
    bRounded = 0;
end

if(~isfield(test, 'note') || ~isfield(test, 'octave'))
    error('Inproper structure fields!');
end

a = 2^(1/12);

reference = [];
reference.note      = 'A';
reference.octave    = 4;
freq_ref            = 440; % Hz

if ischar(test.octave)
    error('test.octave has to be a double')
end

ref_midi    = note2num(reference.note) + reference.octave*12;
test_midi   = note2num(test.note)      + test.octave*12;

halfsteps   = test_midi - ref_midi;
freq        = float('double');
freq        = freq_ref*a^halfsteps;

if bRounded == 1
    freq = round(freq);
end
% fprintf(1, 'Note freq: %.6f\n', freq);
