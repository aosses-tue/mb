function p = setupStandardParamsF0Resynth(p)

p = setupStandardParamsF0mod(p);
p = Ensure_field(p, 'resynth_ftype', 'FIR');
p = Ensure_field(p, 'resynth_ford', 200);
p = Ensure_field(p, 'resynth_order', 'preFilt');
p = Ensure_field(p, 'ramp_atime', 5);
p = Ensure_field(p, 'ramp_rtime', 5);

