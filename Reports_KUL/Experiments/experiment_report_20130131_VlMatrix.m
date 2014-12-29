function experiment_report_20130131_VlMatrix(options)
% function experiment_report_20130131_VlMatrix(options)
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.result_folder   = getpaths('result_folder');

options = Ensure_field(options,'FontSize',16); % Added on 25/08/2014
options = Ensure_field(options, 'bSave', 0);

TrainTrajectory = nan(3,5); % 3 subjects, up to 5 double lists used

% S12, Maria Brughmans: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject_folder          = 'ci-Maria_Brughmans';
stSubject.Initials{1}   = KeepCapitalLetters(subject_folder); % MB
stSubject.filter{1}     = '*lijst*F0m-*';
stSubject.SXX{1}        = 'S12';
stSubject.SNR{1}        = 10;
directory_VlMatrix{1}   = [options.result_folder    subject_folder '/20131125-MT/'];

% S14, Jan Leys: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject_folder          = 'ci-Jan_Leys';
stSubject.Initials{2}   = KeepCapitalLetters(subject_folder); % JL
stSubject.filter{2}     = 'dubbel*F0m-*';
stSubject.SXX{2}        = 'S14';
stSubject.SNR{2}        = 5;
directory_VlMatrix{2}   = [options.result_folder    subject_folder '/20131122-MT/'];

% S17, Julia Schoolmeesters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject_folder          = 'ci-Julia_Schoolmeesters';
stSubject.Initials{3}   = KeepCapitalLetters(subject_folder); % JL
stSubject.filter{3}     = '*lijst*F0m-*';
stSubject.SXX{3}        = 'S17';
stSubject.SNR{3}        = 10;
directory_VlMatrix{3}   = [options.result_folder    subject_folder '/20131126-MT/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = 1:length(stSubject.Initials)

    F0mod_results = dir([directory_VlMatrix{m} stSubject.filter{m} stSubject.Initials{m} '.apr']);

    resF0m.SNR = [];
    for i = 1:length(F0mod_results)
        [xx xx info] = parseapex_VlMatrix([directory_VlMatrix{m} F0mod_results(i).name]);
        
        if str2num(info.SNR) > 50
            info.SNR = 99;  
        else
            info.SNR = str2num(info.SNR);
        end
        
        resF0m.wscores(i)= info.score_words;        
        resF0m.wscores1(i) = info.reliability(1);
        resF0m.wscores2(i) = info.reliability(2);
        resF0m.SNR(i)    = info.SNR;
    end

    idx = find(resF0m.SNR == stSubject.SNR{m});
    
    TrainTrajectory(m,1:length(idx)) = resF0m.wscores(idx);
    
end

clGray = [0.5 0.5 0.5];

figure;
plot(1:length(TrainTrajectory), TrainTrajectory(1,:), 'ko--','MarkerFaceColor','k','LineWidth',2), hold on, grid on
    
plot(1:length(TrainTrajectory), TrainTrajectory(2,:), 's--','Color',clGray,'MarkerFaceColor',clGray,'LineWidth',2)
    
plot(1:length(TrainTrajectory), TrainTrajectory(3,:), 'k>-.','MarkerFaceColor','w','LineWidth',2)
    
Title('F0mod Training trajectory Flemish Matrix')
hl = legend('S12','S14','S17')
set( hl,'FontSize'   , options.FontSize )

ha  = gca;

set(ha,'XTick'      , [1:5] )
set(ha,'FontSize'   , options.FontSize )

Xlabel('Presentation order',options.FontSize)
Ylabel('Word score [%]'    ,options.FontSize)
ylim([30 100])
xlim([0 length(TrainTrajectory)+1])

% set(gcf,'Position',[0 0 1024 300]);
set(gcf, 'PaperPositionMode','auto')

if options.bSave == 1
    saveas(gcf,[options.dest_folder, 'LearningMT.eps'],'epsc');
end
%    