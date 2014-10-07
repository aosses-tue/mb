function [aceData, f0mData, subjects, stLB, test_regs] = figLBScores_xPC(dirExpRes, instr, bPlot, bLoadData, subject_type, Freqs)
% function [aceData, f0mData, subjects, stLB, test_regs] = figLBScores_xPC(dirExpRes, instr, bPlot, bLoadData, subject_type, Freqs)
% 
% stLB  is a structs wich contains all the results of the Loudness Balancing
%       experiments present in the computer
%
% bPlot     = 1; means that the results are going to be plotted
% bLoadData = 1; read the data stored at your computer. If bLoadData = 0,
%       all the loudness balance will be done assuming null values (0 dB
%       correction)
% subject_type = 'NH' or 'CI'; the test_registers are defined depending on
%       the subject_tye.
% 
% Adapted from: figPRScores
%
% Dependencies:
%   readLBData.m
%   freq2note.m
%
% Examples:
%   % Example 1, plotting global results and also the results per subject:
%       close all, 
%       clc, 
%       bPlot = 1; 
%       [deltaACE, deltaF0mod] = figLBScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'UW'}, bPlot);
% 
%   % Example 2, just global result:
%       close all, 
%       clc, 
%       bPlot = 0; 
%       [deltaACE, deltaF0mod] = figLBScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'UW'}, bPlot);
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Starts: ' mfilename])

if nargin < 6
    Freqs = [];
end

if nargin < 5
    subject_type = 'NH';
end

if nargin < 4
    bLoadData = 1;
    disp([mfilename '.m: if you want to avoid data loading, set bLoadData to 0'])
end

if nargin < 3
    bPlot = 0;
end

currentDirectory = cd;
List_of_files = {   'freq2note.m', ...
                    'setupPlotConf.m'};

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stLB                = [];
if strcmp(subject_type,'NH')
    test_regs       = [104, 147, 208, 294];
end
if strcmp(subject_type,'CI')
    test_regs       = [104, 147, 208, 294, 415, 494, 831];
end

if ~isfield(Freqs, 'Fmin')
    Freqs.Fmin      = min(test_regs);
end

if ~isfield(Freqs, 'Fmax')
    stFmax          = freq2note(max(test_regs)); % 4 semitones above Fmax
    Fmax            = sumSemitones2note( stFmax, 4); % 4 semitones above Fmax
    Freqs.Fmax      = round(note2freq(Fmax));
end

k = 1;
Freqs.F(k) = Freqs.Fmin;
stFmin = freq2note(Freqs.Fmin);
newFreq = Freqs.Fmin; % Initial value

while newFreq < Freqs.Fmax
    [xx xx newFreq] = sumSemitones2note(stFmin,k);
    k = k + 1;
    Freqs.F(k) = newFreq;
end

instr = 'UW'; % Experiments will be conducted only using the UW wav files

bFailedReading = 0;
if bLoadData == 1
    try
        [aceData  , f0mData, subjects, ...
         value_ACE, value_F0m] = readLBData(dirExpRes, test_regs, instr);
        if length(test_regs) ~= size(aceData,2)
            bFailedReading = 1;
        end
        if length(test_regs) ~= size(f0mData,2)
            bFailedReading = 1;
            clear aceData;
            clear f0mData;
            clear subjects;
        end
    catch
        bFailedReading = 1;
    end
end

if bFailedReading == 1
    disp([mfilename '.m: data load failed (or dimensions don''t match), check if you have ' subject_type ' data in your computer'])
    disp([mfilename '.m: null amplitudes will be loaded for Loudness Balancing'])
end

if bLoadData == 0 || bFailedReading == 1
    aceData = zeros(1, length(test_regs));
    f0mData = zeros(1, length(test_regs));
    subjects = {'none','none'};
end

for i = 1:size(aceData,1)
    
    Freq_ref        = 131;
    count_less_than_ref = find(test_regs < Freq_ref);
    test_regs_ACE   = [test_regs(1:count_less_than_ref)  Freq_ref   test_regs(count_less_than_ref+1:end)];
    aceData_ACE     = [aceData(i,1:count_less_than_ref)     0       aceData(i,count_less_than_ref+1:end)];
    
    % trick for constant extrapolation at the end: %%%%%%%%%%%%%%%%%%%%%%%%
    aceData_ACE     = [aceData_ACE      aceData_ACE(end)];   
    test_regs_ACE_mod = [test_regs_ACE    max(test_regs_ACE)+10]; 
    f0mData         = [f0mData          f0mData(:,end)];  
    test_regs_mod   = [test_regs        max(test_regs)+10];
    countAdded = size(f0mData, 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display(['Inside ' mfilename ' Caution, reference tone of 131 Hz is being used for interpolation'])
    
    LB_ACE(i,:) = interp1(test_regs_ACE_mod, aceData_ACE ,Freqs.F,'linear','extrap');
    LB_F0m(i,:) = interp1(test_regs_mod    , f0mData(i,:),Freqs.F,'linear','extrap');
    
    f0mData(:,countAdded) = [];
    
%     if bPlot == 1
%         if size(value_ACE,2)<=4
%             figure
%             
%             plot(   Freqs.F, LB_ACE(i,:), 'r-'), hold on
%             plot(   Freqs.F, LB_F0m(i,:), 'b-')
% 
%             legend('ACE', 'F0m'), grid on
%             plot(   test_regs, value_ACE, 'LineWidth', 2, 'Color', 'r', 'Marker', 'o', 'LineStyle','o'), hold on
%             % plot(   test_regs, aceData(i,:), 'LineWidth', 2, 'Color', 'r', 'Marker', 'o', 'LineStyle','o'), hold on
%             plot(   test_regs, value_F0m, 'LineWidth', 2, 'Color', 'b', 'Marker', 'x', 'LineStyle','x')
%             % plot(   test_regs, f0mData(i,:), 'LineWidth', 2, 'Color', 'b', 'Marker', 'x', 'LineStyle','x')
%             title(['Subject: ' name2figname(subjects{i,1}) ' / Loudness Balance'])
%             ylabel('\Delta dB, ref 131 Hz @ 60 dB(A)')
%             xlabel('Frequency [Hz]')
%             ylim([-25 25])
%         end
%     end
    
end

stLB.F = Freqs.F;
stLB.LB_ACE = LB_ACE;
stLB.LB_F0m = LB_F0m;

xTick = test_regs;
yTick = -20:5:10;

if bPlot == 1
    for i=size(aceData,1):-1:1
        idxACE = isnan(aceData(i,1)); % Assumes (for simplicity) that if first value on the row is NaN, the complete row it is NaN
        idxF0m = isnan(f0mData(i,1));
        if idxACE == 1 && idxACE == 1
            aceData(i,:) = [];
            f0mData(i,:) = [];
        end
    end
    
    stPlot.TitleSuffix = {['Loudness Balancing (for ' num2str(size(aceData,1)) ' subjects)']};
    stPlot.xTick = [1:length(xTick)];
    stPlot.xTickLabel = {'104' '147' '208' '294'};
    stPlot.yTick    = yTick;
    stPlot.xLim     = [0 max(stPlot.xTick)+1];
    stPlot.yLim     = [min(yTick)-5 max(yTick)+5];
    stPlot.XLabel   = ['Frequency [Hz]'];
    stPlot.YLabel   = ['\Delta dB ref. UW-131-Hz, ACE'];
    stPlot          = Ensure_field(stPlot,'SeriesLabel',{'xPC F0mod','xPC ACE'});
    stPlot.sameYAxis = 0; % To fix this
    stPlot.LocationLegend = 'NorthEast';
    stPlot          = Ensure_field(stPlot,'bPlotIndividual',0);
    PlotMeans(stPlot, f0mData, aceData);
    
end
cd(currentDirectory);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end

disp(['Ends: ' mfilename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end