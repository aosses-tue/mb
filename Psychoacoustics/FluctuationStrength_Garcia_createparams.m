function FluctuationStrength_Garcia_createparams(N,fs,Ndataset)
% function FluctuationStrength_Garcia_createparams(N,fs,Ndataset)
% 
% Creates the file 'params.mat' (if the file already exists it is deleted
% first) containing all the required parameters for the fluctuation
% strength model.
% 
%   - gzi was deleted on 07/01/2015
%   - Ndataset = 0; is the final fitting, as presented in the thesis
% 
% Author: Rodrigo Garcia
% Original name: Create_params (renamed when copied from Rodrigo's repository)
% Modified by: Alejandro Osses V.
% Last modified: 07/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    N   = 264600;
end
    
if nargin < 3
    Ndataset = 0;
end

OUTPUT_FILE = sprintf('params-%.0f-at-%.0f-Hz-dataset-%.0f.mat',N,fs,Ndataset);
    
% Common parameters:
N0      = round(20 * N / fs) + 1;
N01     = N0 - 1;
Ntop    = round(20e3 * N / fs) + 1;
qb      = N0:1:Ntop;
freqs   = (qb + 1) * fs / N;
Chno    = 47;
zi      = 0.5:0.5:23.5; % Added by AO
Bark        = il_get_Bark;
Barkno      = il_calculate_Barkno(N,fs,qb,Bark);
a0          = il_calculate_a0(N,qb,Barkno);
MinExcdB    = il_calculate_MinExcdB(N0,N01,Ntop,qb,Barkno);
[MinBf zb]  = il_calculate_MinBf(N,N0,fs,Bark,MinExcdB);
Hweight     = il_create_Hweight(N,fs);

switch Ndataset
    case 0
        Cal = 0.1017; % Modified respect to Garcia2015, Eq 6.14
        p_g = 0; % not accounted for in final expression
        p_m = 0.25; % Garcia2015, Eq 6.11
        p_k = 0.375; % Garcia2015, Eq 6.11
        
    case 1
        Cal = 0.1017; % 0.25;
        p_g = 2;
        p_m = 2; 
        p_k = 2;
end
    
if exist(OUTPUT_FILE,'file')
    delete(OUTPUT_FILE);
end

save(OUTPUT_FILE,'-regexp','^(?!OUTPUT_FILE$).*');
 
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bark = il_get_Bark
    Bark = [
        0   0       50      0.5
        1   100     150     1.5
        2   200     250     2.5
        3   300     350     3.5
        4   400     450     4.5
        5   510     570     5.5
        6   630     700     6.5
        7   770     840     7.5
        8   920     1000	8.5
        9   1080	1170	9.5
        10  1270	1370	10.5
        11  1480	1600	11.5
        12  1720	1850	12.5
        13  2000	2150	13.5
        14  2320	2500	14.5
        15  2700	2900	15.5
        16  3150	3400	16.5
        17  3700	4000	17.5
        18  4400	4800	18.5
        19  5300	5800	19.5
        20  6400	7000	20.5
        21  7700	8500	21.5
        22  9500	10500	22.5
        23  12000	13500	23.5
        24  15500   20000   24.5 ]; 
% end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Barkno = il_calculate_Barkno(N,fs,qb,Bark)
    N2      = N / 2 + 1;
    dFs     = fs / N;
    Bark2	= [
        sort([Bark(:,2);Bark(:,3)]),...
        sort([Bark(:,1);Bark(:,4)])
    ];

    Barkno      = zeros(1,N2);
    Barkno(qb)  = interp1(Bark2(:,1),Bark2(:,2),(qb-1)*dFs);
% end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a0 = il_calculate_a0(N,qb,Barkno)
    a0tab = [
        0       0
        10      0
        12      1.15
        13      2.31
        14      3.85
        15      5.62
        16      6.92
        16.5    7.38
        17      6.92
        18      4.23
        18.5    2.31
        19      0
        20      -1.43
        21		-2.59
        21.5	-3.57
        22		-5.19
        22.5	-7.41
        23		-11.3
        23.5	-20
        24		-40
        25		-130
        26		-999
    ];

    a0 = ones(1,N);

    for n = qb;
        a0(n) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(n)));
    end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinExcdB = il_calculate_MinExcdB(N0,N01,Ntop,qb,Barkno)
    HTres = [
        0		130
        0.01    70
        0.17    60
        0.8     30
        1       25
        1.5     20
        2		15
        3.3     10
        4		8.1
        5		6.3
        6		5
        8		3.5
        10		2.5
        12		1.7
        13.3	0
        15		-2.5
        16		-4
        17		-3.7
        18		-1.5
        19		1.4
        20		3.8
        21		5
        22		7.5
        23      15
        24      48
        24.5 	60
        25		130
    ];

    Ntop2   = Ntop - N0 + 1;

    MinExcdB = zeros(1,Ntop2);

    for n =	qb
        MinExcdB(n - N01) = interp1(HTres(:,1),HTres(:,2),Barkno(n));
    end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MinBf zb] = il_calculate_MinBf(N,N0,fs,Bark,MinExcdB)
    Cf = ones(2,24);

    for a = 1:1:24
        Cf(1,a) = round(Bark((a + 1),2) * N / fs) + 1 - N0;
        Cf(2,a) = Bark(a + 1,2);   
    end

    Bf = ones(2,24);
    Bf(1,1) = round(Bark(1,3) * N / fs);
    for a=1:1:24
        Bf(1,a + 1) = round(Bark((a + 1),3) * N / fs) + 1 - N0;
        Bf(2,a)     = Bf(1,a) - 1;
    end
    Bf(2,25) = round(Bark((25),3) * N / fs) + 1 - N0;

    zb      = sort([Bf(1,:),Cf(1,:)]);
    MinBf   = MinExcdB(zb);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hweight = il_create_Hweight(N,fs)

params = struct;
params.sf1 = 0.5;
params.pf1 = 2;
params.pf2 = 8;
params.sf2 = 32;

try
    Hweight = designfilt(...
            'bandpassiir', ...
            'StopbandFrequency1', params.sf1, ...
            'PassbandFrequency1', params.pf1, ...
            'PassbandFrequency2', params.pf2, ...
            'StopbandFrequency2', params.sf2, ...
            'StopbandAttenuation1', 100, ...
            'PassbandRipple', 3, ...
            'StopbandAttenuation2', 100, ...
            'SampleRate', fs);
catch
    d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',params.sf1,params.pf1,params.pf2,params.sf2,80,3,80,fs);
    Hd = design(d,'butter');
    % fvtool(Hd);
    % measure(Hd)

    x = [1; zeros(N-1,1)];
    y = filter(Hd,x);
    freq = (0:(2*pi)/length(x):pi)/pi*fs/2;
    xdft = fft(x);
    ydft = fft(y);

    Hweight = abs(ydft);
    HweightdB = To_dB(Hweight); 
end

% figure;
% plot(freq,20*log10(abs(ydft(1:length(x)/2+1))),'r','linewidth',2);
% xlim([0 50])
% disp('')

% legend('Original Signal','Bandpass Signal');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
