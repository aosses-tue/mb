function new_name = Replace_character(name,Char1,Char2)
% function new_name = Replace_character(name,Char1,Char2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_name        = name;
cont            = find(name == Char1);

if length(cont) ~= 0
    if length(Char1) ~= length(Char2)
        
        Memorise    = name(cont+1:end);
        name(cont:end) = [];
        new_name    = [name Char2 Memorise];
    
    else % Then Char1 and Char2 have the same size
        
        for i = 1:length(cont)
            new_name(cont(i)) = Char2;
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end