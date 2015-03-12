function vocoder_test
% Copyright Tom Francart, 2009

targetpath='/home/tom/temp/audis-cisim/';

[d,fs]=wavread('/mnt/l/speechmaterials/english/bkb/bkbe0101.wav');
r=vocoder(d(:,1),fs);

soundsc(r,fs);
wavwrite(r/max(abs(r)), fs, [targetpath 'bkb101-sim.wav']);
