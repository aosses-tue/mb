function Process_NH_sentence(filename)
% function Process_NH_sentence(filename)
%
% F0mod validation for NH subjects
%
% filename - corresponds to the filename of the Sentence-in-noise results
%
% Strategies:
%   0 = xPC ACE
%   1 = xPC F0mod, IIR Envelope detector
%   2 = xPC F0mod, FIR Envelope detector
%
% % Examples:
%       Process_NH_sentence('NH-Sentence-recognition-retest.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDoSentences    = 1;
bDoWords        = 0;

h = [];
scores_ordered = [];
result_folder   = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';

if nargin == 0
    filename = 'NH-Sentence-recognition.txt';
end

% Sentences:
if bDoSentences == 1
    filename        = [result_folder filename];
    import_scores_from_txt(filename)
    scores          = evalin('base','data');
    scores_txt_info = evalin('base','textdata');

    col_strat   = 1;
    col_scores  = 2;

    idx = find(scores(:,col_strat)==0);
    tmp_strat = scores(idx,col_scores);

    scores_ordered = [scores_ordered tmp_strat(:)];

    idx = find(scores(:,col_strat)==1);
    tmp_strat = scores(idx,col_scores);
    scores_ordered = [scores_ordered tmp_strat(:)];
    
    try
        idx = find(scores(:,col_strat)==2);
        tmp_strat = scores(idx,col_scores);
        scores_ordered = [scores_ordered tmp_strat(:)];
    end

    [p,tbl,stats]           = anova1(scores_ordered);
    [Comparison,Means, hFig]= multcompare(stats);
    h = [h hFig];

    figure(h(end))
    title('Sentence mean analysis')
end

% Words:
if bDoWords == 1
    filename = [result_folder 'NH-Word-recognition.txt'];
    import_scores_from_txt(filename)
    scores_word          = evalin('base','data');
    scores_word_txt_info = evalin('base','textdata');

    col_strat   = 2;
    col_scores  = 3; % wordscores

    idx = find(scores_word(:,col_strat)==0);
    tmp_strat = scores_word(idx,col_scores);

    scores_word_ordered = [];
    scores_word_ordered = [scores_word_ordered tmp_strat(:)];

    idx = find(scores_word(:,col_strat)==1);
    tmp_strat = scores_word(idx,col_scores);
    scores_word_ordered = [scores_word_ordered tmp_strat(:)];

    [p_word,tbl_word,stats_word]           = anova1(scores_word_ordered);
    [Comparison_word,Means_word, hFig]= multcompare(stats_word);
    h = [h hFig];

    figure(h(end))
    title('Word mean analysis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end