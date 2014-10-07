function meas_20130814_implant_type

List_of_files = {   'ConvertSequence_nsb2NMT.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P1 - Sensitivity 12
% P2 - Sensitivity 20
% P3 - Sensitivity 12 + ADRO
% P4 - Sensitivity 12 + ADRO + AGC

filespath = '/home/alejandro/Documenten/Meas/Meas/Calibration-Checking-Map/';
filenames = {   'RP-CIC4-SP12-P1-caltone.out', ...
                'RP-CIC4-SP15-P1-caltone.out', ...
                'RP-CIC3-SP15-P1-caltone.out', ...
                'RP-CIC4-SP12-P2-caltone.out', ...
                'RP-CIC4-SP15-P2-caltone.out', ...
                'RP-CIC3-SP15-P2-caltone.out', ...
                'RP-CIC4-SP12-P3-caltone.out', ...
                'RP-CIC4-SP12-P4-caltone.out'}; 

try
    load( [filespath 'tg_p-q-structures_F0m_Simulink-var'] )
    disp([mfilename '.m: Subject map loaded successfully'])
end

rowPlots = 3;
colPlots = 4;

idxCorrelativeMatrix = [             1:  colPlots ; ...
                            colPlots+1:2*colPlots ; ...
                          2*colPlots+1:3*colPlots ];

limit_y = [0 24];
limit_x = [0 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 1:3;
idxCol = 1;

seq1 = ConvertSequence_txt2nsb([filespath filenames{idx(1)}]);
seq1 = ConvertSequence_nsb2NMT(seq1);

seq2 = ConvertSequence_txt2nsb([filespath filenames{idx(2)}]);
seq2 = ConvertSequence_nsb2NMT(seq2);

seq3 = ConvertSequence_txt2nsb([filespath filenames{idx(3)}]);
seq3 = ConvertSequence_nsb2NMT(seq3);

pulses2plot = min([length(seq1.electrodes) length(seq2.electrodes) length(seq3.electrodes)])

seq1.electrodes         = seq1.electrodes(1:pulses2plot);
seq2.electrodes         = seq2.electrodes(1:pulses2plot);
seq3.electrodes         = seq3.electrodes(1:pulses2plot);
seq1.current_levels     = seq1.current_levels(1:pulses2plot);
seq2.current_levels     = seq2.current_levels(1:pulses2plot);
seq3.current_levels     = seq3.current_levels(1:pulses2plot);

count = find(seq2.electrodes == 24);    seq2.electrodes(count)      = 0;
                                        seq2.current_levels(count)  = 0;
count = find(seq2.electrodes == 24);    seq3.electrodes(count)      = 0;
                                        seq3.current_levels(count)  = 0;

if exist('p','var')
    seq1Patient = MapElectrodogram2PatientBoudaries(seq1,p);
    seq2Patient = MapElectrodogram2PatientBoudaries(seq2,p);
    seq3Patient = MapElectrodogram2PatientBoudaries(seq3,p);
else
    seq1Patient = seq1;
    seq2Patient = seq2;
    seq3Patient = seq3;
end

[y(:,1) hist_y(:,1)] = Average_per_Channel_per_segment(seq1);
[y(:,2) hist_y(:,2)] = Average_per_Channel_per_segment(seq2);
[y(:,3) hist_y(:,3)] = Average_per_Channel_per_segment(seq3);

times2count = 10;
count = hist_y(1:end-1,:) <= 10;
y(count) = 0;

subplot(rowPlots,colPlots,idxCorrelativeMatrix(1,idxCol))
Plot_sequence_for_GUI(seq1Patient);
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(1)} )

subplot(rowPlots,colPlots,idxCorrelativeMatrix(2,idxCol))
Plot_sequence_for_GUI(seq2Patient)
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(2)} )

subplot(rowPlots,colPlots,idxCorrelativeMatrix(3,idxCol))
Plot_sequence_for_GUI(seq3Patient)
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(3)} )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 4:6;
idxCol = 2;

seq1 = ConvertSequence_txt2nsb([filespath filenames{idx(1)}]);
seq1 = ConvertSequence_nsb2NMT(seq1);

seq2 = ConvertSequence_txt2nsb([filespath filenames{idx(2)}]);
seq2 = ConvertSequence_nsb2NMT(seq2);

seq3 = ConvertSequence_txt2nsb([filespath filenames{idx(3)}]);
seq3 = ConvertSequence_nsb2NMT(seq3);

pulses2plot = min([length(seq1.electrodes) length(seq2.electrodes) length(seq3.electrodes)])

seq1.electrodes         = seq1.electrodes(1:pulses2plot);
seq2.electrodes         = seq2.electrodes(1:pulses2plot);
seq3.electrodes         = seq3.electrodes(1:pulses2plot);
seq1.current_levels     = seq1.current_levels(1:pulses2plot);
seq2.current_levels     = seq2.current_levels(1:pulses2plot);
seq3.current_levels     = seq3.current_levels(1:pulses2plot);

count = find(seq2.electrodes == 24);    seq2.electrodes(count)      = 0;
                                        seq2.current_levels(count)  = 0;
count = find(seq2.electrodes == 24);    seq3.electrodes(count)      = 0;
                                        seq3.current_levels(count)  = 0;

if exist('p','var')
    seq1Patient = MapElectrodogram2PatientBoudaries(seq1,p);
    seq2Patient = MapElectrodogram2PatientBoudaries(seq2,p);
    seq3Patient = MapElectrodogram2PatientBoudaries(seq3,p);
else
    seq1Patient = seq1;
    seq2Patient = seq2;
    seq3Patient = seq3;
end

[y(:,1) hist_y(:,1)] = Average_per_Channel_per_segment(seq1);
[y(:,2) hist_y(:,2)] = Average_per_Channel_per_segment(seq2);
[y(:,3) hist_y(:,3)] = Average_per_Channel_per_segment(seq3);

times2count = 10;
count = hist_y(1:end-1,:) <= 10;
y(count) = 0;

subplot(rowPlots,colPlots,idxCorrelativeMatrix(1,idxCol))
Plot_sequence_for_GUI(seq1Patient);
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(1)} )

subplot(rowPlots,colPlots,idxCorrelativeMatrix(2,idxCol))
Plot_sequence_for_GUI(seq2Patient)
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(2)} )

subplot(rowPlots,colPlots,idxCorrelativeMatrix(3,idxCol))
Plot_sequence_for_GUI(seq3Patient)
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(3)} )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 7:8;
idxCol = 3;

seq1 = ConvertSequence_txt2nsb([filespath filenames{idx(1)}]);
seq1 = ConvertSequence_nsb2NMT(seq1);

seq2 = ConvertSequence_txt2nsb([filespath filenames{idx(2)}]);
seq2 = ConvertSequence_nsb2NMT(seq2);

pulses2plot = min([length(seq1.electrodes) length(seq2.electrodes) ])

seq1.electrodes         = seq1.electrodes(1:pulses2plot);
seq2.electrodes         = seq2.electrodes(1:pulses2plot);
seq1.current_levels     = seq1.current_levels(1:pulses2plot);
seq2.current_levels     = seq2.current_levels(1:pulses2plot);

count = find(seq2.electrodes == 24);    seq2.electrodes(count)      = 0;
                                        seq2.current_levels(count)  = 0;


if exist('p','var')
    seq1Patient = MapElectrodogram2PatientBoudaries(seq1,p);
    seq2Patient = MapElectrodogram2PatientBoudaries(seq2,p);
else
    seq1Patient = seq1;
    seq2Patient = seq2;
end

[y(:,1) hist_y(:,1)] = Average_per_Channel_per_segment(seq1);
[y(:,2) hist_y(:,2)] = Average_per_Channel_per_segment(seq2);

times2count = 10;
count = hist_y(1:end-1,:) <= 10;
y(count) = 0;

subplot(rowPlots,colPlots,idxCorrelativeMatrix(1,idxCol))
Plot_sequence_for_GUI(seq1Patient);
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(1)} )

subplot(rowPlots,colPlots,idxCorrelativeMatrix(2,idxCol))
Plot_sequence_for_GUI(seq2Patient)
h = gca;
set(h,'YTick',1:2:23); set(h,'YTickLabel',22:-2:0)
axis([limit_x limit_y])
title( filenames{idx(2)} )

end