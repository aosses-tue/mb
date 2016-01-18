function Play_random_one_ref_two_test_tones(file_ref,file_test1, file_test2, Ntimes, dur)
% function Play_random_one_ref_two_test_tones(file_ref,file_test1, file_test2, Ntimes, dur)
%
% 1. Description:
%
% 2. Stand-alone example:
%       dir = 'C:\Users\aosses\Desktop\group4\96000-Hz-24-bits-all-files\';
%       fref = [dir '24-bit-flac-96000-Hz-24-b.wav'];
%       ftest1 = [dir '128-mp3-96000-Hz-24-b.wav'];
%       ftest2 = [dir '320-mp3-96000-Hz-24-b.wav'];
%       dur = 5;
%       Ntimes = 10;
%       Play_random_one_ref_two_test_tones(fref,ftest1, ftest2, Ntimes, dur);
% 
%       condition = 1;
%       dir = ['C:\Users\aosses\Desktop\group4\96000-Hz-24-bits-filtered\LP-cond-' num2str(condition) delim];
%       fref = [dir '24-bit-flac-96000-Hz-24-b.wav'];
%       ftest1 = [dir '128-mp3-96000-Hz-24-b.wav'];
%       ftest2 = [dir '320-mp3-96000-Hz-24-b.wav'];
%       dur = 5;
%       Ntimes = 10;
%       Play_random_one_ref_two_test_tones(fref,ftest1, ftest2, Ntimes, dur);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 05/01/2016
% Last update on: 05/01/2016 
% Last use on   : 05/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xref fs ] = Wavread(file_ref);
[xt1  fs1] = Wavread(file_test1);
[xt2  fs2] = Wavread(file_test2);

if fs ~= fs1 | fs ~= fs2
    error('Audio files should have the same sample frequency');
end

if nargin == 5
    idxdur = round(dur*fs);
    
    if size(xref,1) > idxdur
        xref = xref(1:idxdur,:);
    end
    
    if size(xt1,1) > idxdur
        xt1 = xt1(1:idxdur,:);
    end
    
    if size(xt2,1) > idxdur
        xt2 = xt2(1:idxdur,:);
    end
end

durs(1) = size(xref,1)/fs;
durs(2) = size(xt1,1)/fs;
durs(3) = size(xt2,1)/fs;

for i = 1:Ntimes
    
    bRefFirst   = randi(2) - 1;
    bTest1first = randi(2) - 1;
    if bRefFirst
        sound(xref,fs);
        pause(durs(1)*1.05);

        label1Sound = file_ref;

        if bTest1first
            sound(xt1,fs);
            label2Sound = file_test1;
        else
            sound(xt2,fs);
            label2Sound = file_test2;
        end
    else
        if bTest1first
            sound(xt1,fs);
            label1Sound = file_test1;
            pause(durs(2)*1.05)
        else
            sound(xt2,fs);
            label1Sound = file_test2;
            pause(durs(3)*1.05)
        end
        sound(xref,fs);
        label2Sound = file_ref;
    end

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Give your answer (write it down somewhere) and then press any button')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause()

    disp('The correct answer was:');
    fprintf('Attempt No %.0f:\n\tFirst sound was: %s\n\tSecond sound was: %s\n',i,label1Sound,label2Sound);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
