function [p, q, Max_corr] = meas_20121205
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [p, q, Max_corr] = meas_20121205.m
%
% Synchronises 2 electrodograms measuring the same experiment and being
% monitored using the DIET CIC4 Decoder.
% These measurements (made on Dec 5th, 2012) are not meaningful...just
% first approaches...
%
% Load the following files (measured with DIET CIC4):
%       test_XPC_05.out -> Previously saved as test_XPC_05.mat
%       test_XPC_06.out -> Previously saved as test_XPC_06.mat
%       test_XPC_10.out -> Previously saved as test_XPC_10.mat
%
% Dependencies:
%       NIC_analyse_data.m
%       Subplot_sequence.m
%       test_XPC_05.mat
%       test_XPC_06.mat
%       test_XPC_10.mat
%       PureTones_1_2_4kHz_at_16kHz.wav
%       PureTones_1_2_4kHz_at_16kHz+Noise.wav
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Programmed by AOV, ExpORL, KULeuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumMaxima = 8;

display('Start script "Electrodogram_Sync.m"')
display('Part 1: Synchronisation of test_XPC_05 and _06')

load('test_XPC_05.mat');
load('test_XPC_06.mat');
load('test_XPC_10.mat');

y05 = NIC_analyse_data(test_XPC_05, NumMaxima); 
y06 = NIC_analyse_data(test_XPC_06, NumMaxima); 
y10 = NIC_analyse_data(test_XPC_10, NumMaxima); 

Figure2 = figure; 
cant_subplot = 3; 
Subplot_sequence(y05, 'Simulink', [22+1-(0:22)],0,Figure2, 1,'v',cant_subplot), title('05')
Subplot_sequence(y06, 'Simulink', [22+1-(0:22)],0,Figure2, 2,'v',cant_subplot), title('06')
Subplot_sequence(y10, 'Simulink', [22+1-(0:22)],0,Figure2, 3,'v',cant_subplot), title('10')

% y10 is saturated, I will use for the correlation tests signals 05 and 06

q05 = y05;
q06 = y06;

y05.current_levels_matrix = buffer( y05.current_levels, NumMaxima);
y05.electrodes_matrix = buffer( y05.electrodes, NumMaxima);

y06.current_levels_matrix = buffer( y06.current_levels, NumMaxima);
y10.current_levels_matrix = buffer( y10.current_levels, NumMaxima);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[size_row size_col] = size( y05.current_levels_matrix );

CompChannel = 15; % Channel to be correlated

y05.RowTest     = zeros( 1, size_col );
cont            = find(y05.electrodes==CompChannel);
cont_int        = floor(cont/8)+1;
y05.RowTest( cont_int ) = y05.current_levels( cont ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[size_row size_col] = size( y06.current_levels_matrix );

y06.RowTest     = zeros( 1, size_col );
cont            = find(y06.electrodes==CompChannel);
cont_int        = floor(cont/8)+1;
y06.RowTest( cont_int ) = y06.current_levels( cont ); 

%%%% Correlation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaPulsChann = 1/900; % Number of pulses per sec (7200 pps in 8 channels)

if length(y05.RowTest) ~= length(y06.RowTest)
    y0_corr = (xcorr(y05.RowTest, y06.RowTest,'none'));
else
    y0_corr = (xcorr(y05.RowTest, y06.RowTest,'coeff'));
end

[max_a max_b] = max(y0_corr);

Delay = (length(y0_corr) + 1)/2 - max_b;

if Delay > 0
    display('y05 then y06')
else
    display('y06 then y05')
end

Delay*deltaPulsChann % Then I have to delay the electrodes in 8*Delay
y05sync = q05;
y06sync = q06;
y05sync.electrodes = y05.electrodes(1:19992 - Delay*8);
y05sync.current_levels = y05.current_levels(1:19992 - Delay*8);

y06sync.electrodes = y06.electrodes(1+Delay*8 : 19992);
y06sync.current_levels = y06.current_levels(1:19992 - Delay*8);

Figure2 = figure; 
cant_subplot = 2; 
Subplot_sequence(y05sync, 'Simulink', [22+1-(0:22)],0,Figure2, 1,'v',cant_subplot), title('Sync 05')

Subplot_sequence(y06sync, 'Simulink', [22+1-(0:22)],0,Figure2, 2,'v',cant_subplot), title('Sync 06')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Part 2: Synchronisation of "PureTones_1_2_4kHz_at_16kHz.wav" and *.+Noise.wav')

[x1 fs1] = wavread('PureTones_1_2_4kHz_at_16kHz.wav');
[x2 fs2] = wavread('PureTones_1_2_4kHz_at_16kHz+Noise.wav');

deltaT = 1/fs1; 

if length(x1) ~= length(x2)
    y_corr = (xcorr(x1, x2,'none'));
else
    y_corr = (xcorr(x1, x2,'coeff'));
end

[max_a max_b] = max(y_corr);

retardo = (length(y_corr) + 1)/2 - max_b;

if retardo > 0
    display('x_1 then x_2')
else
    display('x_2 then x_1')
end

retardo*deltaT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('End script "Electrodogram_Sync.m"')