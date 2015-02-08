% Received on 20150206 from AK (e-mail)

clear
clf
inputfile = 't54.wav';
[s1, sf] = wavread(inputfile);
compressionFactor = 4;
t1 = [1:length(s1)]/sf;
s2 = s1(1:compressionFactor:length(s1));
sound(s2,sf/compressionFactor)
wavwrite(s2, sf/compressionFactor, 16, 't54_2.wav')
t2 = [1:length(s2)]/(sf/compressionFactor);
subplot(2,1,1)
plot(t1, s1)
axis([0.1*t1(length(t1)) 0.12*t1(length(t1)) -1 1])
subplot(2,1,2)
plot(t2, s2)
axis([0.1*t2(length(t2)) 0.12*t2(length(t2)) -1 1])

disp('')