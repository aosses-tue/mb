function p = setupStandardParamsF0mod(p)

if nargin == 0
   p = []; 
end
p = setupStandardParams(p);
p = Ensure_field(p, 'f0File', 0);
p = Ensure_field(p, 'plot_modulator', 0);