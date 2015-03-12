
function [output_time output_spectra] = STI_OctFilterBank(x,fs,midfreq)
% input = 'DT1_short.wav';
% 
% [x,fs]=wavread(input);
% D = 4;
% x = x(1:D:end)';
% fs = fs/D;

% x = input;
% x = real(ifft(ones(1,1000)));
% fs = 2000;
X = (fft(x));
N = length(x);
freq = linspace(0,fs/2,N/2 +1);
f_axis = [-1*freq(end:-1:2)  freq(1:end-1)];

%resolution of data
resol=freq(2)-freq(1);

%band cross-over frequencies
crossfreq(1)=midfreq(1)/(2^(1/2));
crossfreq(2:length(midfreq)+1)=midfreq*(2^(1/2));

%cross-over indicies
y=crossfreq/resol;

%initialize output matrix
output_time = zeros(N,length(midfreq));
output_spectra = zeros(N,length(midfreq));
%rounding up
crosselem=round(y);
for n=1:length(y)
if crosselem(n)<y(n)
    crosselem(n)=crosselem(n)+1;
end
end

nn=1;

while crossfreq(nn+1)<=freq(end)
    
    output_spectra(:,nn) = scut(X,crossfreq(nn),crossfreq(nn+1),fs);
    output_time(:,nn) = real(ifft(output_spectra(:,nn)));
    nn=nn+1;
     if 1+nn > length(crossfreq)
        break
    end
end

