function AnalyzeAudioviz

fs=48000;
% part1file='~/dvdrip-data/audioviz/vob/001/rawaudio001.wav';
part1file='~/dvdrip-data/audioviz/vob/001/rawaudio-all.wav';



for q=1:2
    
    if (q==1)       % inge
        totalfrom=0;            % start analyzing from here
        totalto=30*60;
    else
        totalfrom=31*60;
        totalto=60*60;
    end

    % totalsamples=wavread(part1file,'size');
    totalsamples=(totalto-totalfrom)*fs;

    partlen=20;         % seconds
    N=floor(totalsamples(1)/fs/partlen);               % number of parts

    resolution=20;
    % threshold=1e-3;
    threshold=0;

    NFFT = 2^ceil(log(fs/resolution)/log(2));    %The number of points in the FFT
    wind = boxcar(NFFT);  
    Pxx = zeros(NFFT/2+1,N);               
    rmss = zeros(1,N);

    for i=1:N
        disp(['Part ' num2str(i)]);
        from=(totalfrom+(i-1)*partlen)*fs+1;
        to=(totalfrom+i*partlen)*fs;

        [audio,rate] = wavread(part1file, [from to]);
        audio=audio(:,1);
        if (rate~=fs)
            error('Invalid rate');
        end

        %Remove all silence portions out of the wav file
        if threshold > 0
            audio = removegap(audio,rate,threshold);    
        end

        [Pxx(:,i),F] = psd(audio,NFFT,rate,wind,0,'none');  %Calculate the average Power Spectrum Density of the wav file
        Nsections(i) = ceil(length(audio)/NFFT);        %The number of blocks (of NFFT samples) that are in the wavfile
        Pxx(:,i) = Nsections(i)*Pxx(:,i);               %Weigth the importance of this wav files Pxx by Nsections
        
        rmss(i) = meanrms_thresh(audio,fs,0.001);
    end

    Pxx_average = sum(Pxx,2)/sum(Nsections);
    rms_average = mean(rmss(~isnan(rmss)));
    
    if (q==1)
        Pxx_inge=Pxx_average;
        Rms_inge=rms_average;
    else
        Pxx_man=Pxx_average;
        Rms_man=rms_average;
    end


end

Pxx_both=mean( [Pxx_inge Pxx_man], 2 );

figure;
semilogx(F,10*log10(Pxx_inge),F,10*log10(Pxx_man), F, 10*log10(Pxx_both));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Average spectra for the audiovisual sentences');
legend('Spectrum Inge', 'Spectrum man', 'Average');

disp(sprintf('RMS inge=%2.2f', 20*log10(Rms_inge)));
disp(sprintf('RMS man =%2.2f', 20*log10(Rms_man)));
disp(sprintf('RMS avg =%2.2f', 20*log10(mean([Rms_man Rms_inge])));

save AudiovizSpectra;


%Helper function to remove the silenced gaps in the between words
function y = removegap(x,fs,thr);
blocklength = round(0.02*fs);  %Use a 20ms window
VA = zeros(size(x));
for i = 1:length(x)
    VA(i) = (sqrt(mean(x(i:min(i+blocklength-1,length(x))).^2)) > thr);
end
y = x(VA);
return;