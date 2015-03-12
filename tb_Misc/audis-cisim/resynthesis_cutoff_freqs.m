function f=resynthesis_cutoff_freqs(nbands, type)
% Copyright Tom Francart, 2009

switch type
    case 'greenwood'
        insertion_depth=22;     % insertion depth of electrode (mm)
        cochlea_length=33;      % length of cochlea (mm)
        array_length=18;
        
        elec_spacing=array_length/nbands;
        positions=[cochlea_length-insertion_depth:elec_spacing:cochlea_length-insertion_depth+array_length];
        f=Greenwood_x2cf(positions)';
        
    otherwise
        error('Invalid type');
end
