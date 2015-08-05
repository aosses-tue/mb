function p = setupStandardParamsACEResynth(p)
% function p = setupStandardParamsACEResynth(p)
%
% Programmed by Matthias Milczynski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = setupStandardParams(p);
p = Ensure_field(p, 'resynth_ftype', 'FIR');
p = Ensure_field(p, 'resynth_ford', 200);
p = Ensure_field(p, 'resynth_order', 'preFilt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end