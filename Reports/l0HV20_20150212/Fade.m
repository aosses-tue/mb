function env = envelope(s, sf, rampDuration)
env = ones(size(s));
fadein = (0:1/sf: rampDuration)/ rampDuration;
env(1:length(fadein)) = fadein;
env(length(s):-1:length(s)-length(fadein)+1) = fadein;
env = s.*env;
