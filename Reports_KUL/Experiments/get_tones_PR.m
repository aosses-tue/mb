function [tones, freqs] = get_tones_PR(Registers, bIncludeReference)

% function [tones, freqs] = get_tones_PR(Registers, bIncludeReference)
%
% Inputs:
%   Registers
%   bIncludeReference:  if 1, the first column will contain the Reference tone
%
% Programmed by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    bIncludeReference = 1;
end

for i = 1:length(Registers)
    
    Reg = Registers{i};
    count = find(Reg == '_'); 
    Reg(count)  = []; % Deleting '_' in case it is present in Registers
    
    tones{i,1}  = Reg;
    
    Note.octave = Reg(end);
    Note.note   = Reg(1:end-1);
    
    [a, tones{i,1}, freqs(i,1)] = sumSemitones2note(Note,0);
    [a, tones{i,2}, freqs(i,2)] = sumSemitones2note(Note,1);
    [a, tones{i,3}, freqs(i,3)] = sumSemitones2note(Note,2);
    [a, tones{i,4}, freqs(i,4)] = sumSemitones2note(Note,3);
    [a, tones{i,5}, freqs(i,5)] = sumSemitones2note(Note,4);
    
end

if bIncludeReference == 0
    tones(:,1) = []; % We eliminate the reference column
    freqs(:,1) = []; % We eliminate the reference column
end