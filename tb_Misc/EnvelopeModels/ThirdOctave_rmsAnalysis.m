function rms_out = ThirdOctave_rmsAnalysis(x,fs,midfreq)
% function rms_out = ThirdOctave_rmsAnalysis(x,fs,midfreq)
%
% Original name: STI_ThirdOctave_rmsAnalysis
%-----------------------------------------
%AVERAGING IN 1/3-OCTAVE BANDS, 
%-----------------------------------------

N = length(x);
X = (fft(x));
X_mag  = abs(X) ;
X_power = X_mag.^2/N ;% power spectrum.
X_power_pos = X_power(1:fix(N/2)+1) ;
X_power_pos(2:end) = X_power_pos(2:end).* (2)  ; %take positive frequencies only and mulitply by two-squared to get the same total energy(used since the integration is only performed for positive freqiencies)

freq= linspace(0,fs/2,length(X_power_pos));

%resolution of data
resol=freq(2)-freq(1);

%band cross-over frequencies

if nargin<3
    [0.63 .8 1 1.25 1.6 2.0 2.5 3.15 4.0 5.0 6.3 8 10 12.5];
end

crossfreq(1)=midfreq(1)/(2^(1/6));
crossfreq(2:length(midfreq)+1)=midfreq*(2^(1/6));

%cross-over indicies
y=crossfreq/resol;

crosselem=round(y);
for n=1:length(y)
    if crosselem(n)<y(n)
        crosselem(n)=crosselem(n)+1;
    end
end

nn=1;

while crossfreq(nn+1)<=freq(end)% for nn =1:length(crossfreq)-1    
    rms_out(nn) = sqrt(sum(X_power_pos(crosselem(nn):crosselem(nn+1)-1))/N);
    % rms_out(nn) = (max(X_power_pos(crosselem(nn):crosselem(nn+1)-1))/N);
   
    nn=nn+1;  
    if 1+nn > length(crossfreq)
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
