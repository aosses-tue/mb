function output = exp_fastl2007(flag)
%EXP_FASTL2007  Figures from Fastl and Zwicker (2007)
%   Usage: output = exp_fastl2007(flag);
%
%   EXP_FASTL2007 reproduces Fluctuation strength data from the book Fastl 
%       and Zwicker (2007).
%
%   Examples:
%   ---------
%     exp_fastl2007('fmod');
%
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 02/08/2015
% Last update on: 02/08/2015 
% Last use on   : 08/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    flag = 3;
end

where_data = 'D:\MATLAB_shared\Psychoacoustics\FluctuationStrength\Literature\Data-processed\';

switch flag
    case {0,'fmod'}
        
        [FS_AM_std1, FS_AM_std2] = il_scan_file([where_data 'AM-fmod.csv']);
        % % Approximated estimation:
        % FS_AM_std1 = [3   10.25 26.25 66   100  99 30    9.5; ...
        %   0    7.75 20    50    98  75  7    0; ...
        %   5.5 20    50    90.5 102 110 49.5 12.25];

        % % Approximated estimation:
        % FS_AM_std2 = [5.5 15 33 53.5 91.5 99 25 0.25; ...
        %               1 12 26 39.5 65.5 73 14 -3; ...
        %               8.25 16.25 36 65.5 109.5 110 58 12];

        [FS_FM_std1, FS_FM_std2] = il_scan_file([where_data 'FM-fmod.csv']);

        % FS_FM_std1 = [10.2 22.75 35 57.25 100 50   19.5   0.75; ...
        %                5.5  5    21 50     98 30.5 11    -1.5; ...
        %                2   42    79 95    101 73   39.75  9.75];
        % FS_FM_std2 = [10.25 20.5  50.25 60.5 100    20.5   7.5  2; ...
        %                8    19.25 40.5  55.5  75.75 10.25  4    0; ...
        %               13.5  23    60.25 89.5 110    90    35   10];
        
        xtest_AM = [0.25 0.5 1 2 4 8 16 32]; % Hz
        xtest_FM = [0.25 0.5 1 2 4 8 16 32]; % Hz
        xtest_type = 'fmod'; % Hz
        text_XLabel = 'f_m_o_d [Hz]';
        
        FS_BBN_std1 = [22.5  57   79   101  100    82   29     6; ...
                        3.5  29.5 59.5  90   97    68.5 20     1; ... 
                       30    70   98   110  102.25 84.5 44    14];
        FS_BBN_std2 = [12.5  27   45.5  75   99.75 92   24.25  9.75; ...
                        5.5  24.5 39    51.5 60.5  67   13     4; ...
                       16    30   48    92  110   110   32    19.5]; 
        output.FS100 = [1.3 2 1.8]; % vacil, 100% of fluctuation
        
    case {1,'fc'}
        [FS_AM_std1, FS_AM_std2] = il_scan_file([where_data 'AM-fc.csv']);
        [FS_FM_std1, FS_FM_std2] = il_scan_file([where_data 'FM-fc.csv']);
        FS_BBN_std1 = [];
        FS_BBN_std2 = [];
        
        xtest_AM = [125  250  500 1000 2000 4000 8000]; % Hz
        xtest_FM = [500 1000 1500 2000 3000 4000 6000 8000]; % Hz
        
        xtest_type = 'fc'; % Hz
        text_XLabel = 'f_c [Hz]';
        
        output.FS100 = [1.4 2.22 NaN]; % vacil, 100% of fluctuation
        
    case {2,'SPL'}
        [FS_AM_std1, FS_AM_std2] = il_scan_file([where_data 'AM-SPL.csv'],3);
        [FS_FM_std1, FS_FM_std2] = il_scan_file([where_data 'FM-SPL.csv'],3);
        FS_BBN_std1 = [];
        FS_BBN_std2 = [];
        
        xtest_AM = [50 60 70 80 90]; % Hz
        xtest_FM = [40 50 60 70 80]; % Hz
        xtest_type = 'SPL'; % dB
        text_XLabel = 'SPL [dB]';
        
        output.FS100 = [2.4 1.85 3]; % vacil, 100% of fluctuation
        
    case {3,'md','df'}
        xtest_AM = [0  3.2 6.8 9.1 11 13 16 19 40]; % Hz
        xtest_FM = [16 32 100 200 300 500 700]; % Hz
                
        [FS_AM_std1, FS_AM_std2] = il_scan_file([where_data 'AM-md.csv'],3);
        [FS_FM_std1, FS_FM_std2] = il_scan_file([where_data 'FM-df.csv'],3);
        xtest_type = {'Modulation depth [dB]','Frequency deviation [Hz]'}; % dB
        FS_BBN_std1 = [];
        FS_BBN_std2 = [];
        
        output.FS100 = [1.35 2 1.75]; % vacil, 100% of fluctuation
        
end

output.xtest_AM = xtest_AM; 
output.xtest_FM = xtest_FM; 
output.xtest_type = xtest_type;
output.FS_BBN_std1 = FS_BBN_std1;
output.FS_BBN_std2 = FS_BBN_std2;
output.FS_AM_std1  = FS_AM_std1;
output.FS_AM_std2  = FS_AM_std2;
try
    output.FS_AM = prctile([FS_AM_std1(1,:);FS_AM_std2(1,:)],50);
    output.FS_FM = prctile([FS_FM_std1(1,:);FS_FM_std2(1,:)],50);
catch % then standard 2 is not available from data
    output.FS_AM = FS_AM_std1;
    output.FS_FM = FS_FM_std1;
end

output.FS_FM_std1  = FS_FM_std1;
output.FS_FM_std2  = FS_FM_std2;

if nargout == 0
    
    dx = 0.15;
    text_YLabel = 'Relative FS [%]';
    Ntest = length(xtest_AM);
        
    figure; 
    subplot(1,3,1) % AM-BBN
    iq1L = FS_BBN_std1(1,:) - FS_BBN_std1(2,:);
    iq1U = FS_BBN_std1(3,:) - FS_BBN_std1(1,:);
    iq2L = FS_BBN_std2(1,:) - FS_BBN_std2(2,:);
    iq2U = FS_BBN_std2(3,:) - FS_BBN_std2(1,:);
    errorbar([1:Ntest]-dx,FS_BBN_std1(1,:),iq1L,iq1U,'bx','LineWidth',2); hold on
    errorbar([1:Ntest]+dx,FS_BBN_std2(1,:),iq2L,iq2U,'ro');
    plot([1:Ntest],mean([FS_BBN_std1(1,:); FS_BBN_std2(1,:)]),'k>--','LineWidth',2)
    xlabel(text_XLabel)
    ylabel(text_YLabel)
    grid on;
    xlim([0.5 Ntest+0.5])
    ha = gca;
    set(ha,'XTick',[1:Ntest])
    set(ha,'XTickLabel',xtest_AM)
    title('AM BBN')
 
    subplot(1,3,2) % AM-SIN
    iq1L = FS_AM_std1(1,:) - FS_AM_std1(2,:);
    iq1U = FS_AM_std1(3,:) - FS_AM_std1(1,:);
    iq2L = FS_AM_std2(1,:) - FS_AM_std2(2,:);
    iq2U = FS_AM_std2(3,:) - FS_AM_std2(1,:);
    errorbar([1:Ntest]-dx,FS_AM_std1(1,:),iq1L,iq1U,'bx','LineWidth',2); hold on
    errorbar([1:Ntest]+dx,FS_AM_std2(1,:),iq2L,iq2U,'ro');
    plot([1:Ntest],mean([FS_AM_std1(1,:); FS_AM_std2(1,:)]),'k>--','LineWidth',2)
    xlabel(text_XLabel)
    ylabel(text_YLabel)
    grid on;
    xlim([0.5 Ntest+0.5])
    ha(end+1) = gca;
    set(ha(end),'XTick',[1:Ntest])
    set(ha(end),'XTickLabel',xtest_AM)
    title('AM SIN')
    
    subplot(1,3,3) % AM-SIN
    iq1L = FS_FM_std1(1,:) - FS_FM_std1(2,:);
    iq1U = FS_FM_std1(3,:) - FS_FM_std1(1,:);
    iq2L = FS_FM_std2(1,:) - FS_FM_std2(2,:);
    iq2U = FS_FM_std2(3,:) - FS_FM_std2(1,:);
    errorbar([1:Ntest]-dx,FS_FM_std1(1,:),iq1L,iq1U,'bx','LineWidth',2); hold on
    errorbar([1:Ntest]+dx,FS_FM_std2(1,:),iq2L,iq2U,'ro');
    plot([1:Ntest],mean([FS_FM_std1(1,:); FS_FM_std2(1,:)]),'k>--','LineWidth',2)
    xlabel(text_XLabel)
    ylabel(text_YLabel)
    grid on;
    xlim([0.5 Ntest+0.5])
    ha(end+1) = gca;
    set(ha(end),'XTick',[1:Ntest])
    set(ha(end),'XTickLabel',xtest_FM)
    title('FM SIN')
    
    legend('4-Hz std','0.5-Hz std','avg')
    linkaxes(ha,'xy');
    ylim([-5 110])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Std1, Std2] = il_scan_file(filename,cols2read)

if nargin < 2
    cols2read = 5;
end

alltext = Get_text_from_txt(filename);
        
r = textscan(alltext{1},'%s',cols2read,'Delimiter',',');
r = r{1};
j1 = 1;
j2 = 1;
col_l = [];
col_u = [];
Std1 = [];
Std2 = [];

for i = 1:cols2read
    if strcmp( r(i) ,'m')
        col_median = i;
    end

    if strcmp( r(i) ,'l')
        col_l = i;
    end

    if strcmp( r(i) ,'u')
        col_u = i;
    end
end

for i = 2:length(alltext)
    r = textscan(alltext{i},'%s',cols2read,'Delimiter',',');
    r = r{1};
    switch r{1}
        case {'1'}
            Std1(1,j1) = str2num(r{col_median});
            if length(col_l) ~= 0
                Std1(2,j1) = str2num(r{col_l});
                Std1(3,j1) = str2num(r{col_u});
            end
            j1 = j1 + 1;

        case {'2'}
            Std2(1,j2) = str2num(r{col_median});
            if length(col_l) ~= 0
                Std2(2,j2) = str2num(r{col_l});
                Std2(3,j2) = str2num(r{col_u});
            end
            j2 = j2 + 1;
    end
end
