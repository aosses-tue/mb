function y = r20160511_learning_multidimensional_space
% function y = r20160511_learning_multidimensional_space
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 11/05/2016
% Last update on: 11/05/2016 
% Last use on   : 11/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% http://nl.mathworks.com/help/stats/examples/classical-multidimensional-scaling.html

close all

% X(1,1:8) are the treatments
X = [39.1     18.7;
     40.7     21.2;
     41.5     21.5;
     39.2     21.8;
     38.7     20.6;
     41.7     20.1;
     40.1     22.1;
     39.2     21.6];

% D would be the Similarity score
D = [4.69 6.79 3.50 3.11 4.46 5.57 3.00 ...
          2.10 2.27 2.65 2.36 1.99 1.74 ...
               3.78 4.53 2.83 2.44 3.79 ...
                    1.98 4.35 2.07 0.53 ...
                         3.80 3.31 1.47 ...
                              4.35 3.82 ...
                                   2.57];
                               
D = D/max(D);
squareform(D)

% Y is the resulting multidimensional space ('configuration matrix')
%   - size(Y) = n x p
%   - rows of Y: coordinates of n points in the p-dimensional space
[Y,eigvals] = cmdscale(D);

dim2account = 2;
Yred = Y(:,1:dim2account);

% The best fit would not be used at this moment:
[D,Z] = procrustes(X,Yred); % distance and best fit for multidimensional scaling (only 2D)

% figure;
% plot(X(:,1),X(:,2),'bo',Z(:,1),Z(:,2),'rd');
% labels = num2str((1:8)');
% text(X(:,1)+.05,X(:,2),labels,'Color','b');
% text(Z(:,1)+.05,Z(:,2),labels,'Color','r');
% xlabel('Distance East of Reference Point (Km)');
% ylabel('Distance North of Reference Point (Km)');
% legend({'Spatial Locations','Constructed Genetic Locations'},'Location','SE');

xoffset = 0.025;

figure;
plot( Yred(:,1),Yred(:,2),'rd' ); grid on
labels = num2str((1:8)');
% text(X(:,1)+.05,X(:,2),labels,'Color','b');
text(Yred(:,1)+xoffset,Yred(:,2),labels,'Color','r');
xlabel('Dimension 1');
ylabel('Dimension 2');
% legend({'2D space'},'Location','SE');

pairs = Get_pairwise_combinations(1,8);

for i = 1:size(pairs,1)
    data1 = Yred(pairs(i,1),:);
    data2 = Yred(pairs(i,2),:);
    leg4plot{i} = sprintf('%.0f%.0f',pairs(i,1),pairs(i,2));
    distance(i,1) = il_euclidean_dist( data1,data2 );
end

% idx2plot = 1:7;
idx2plot  = find(pairs(:,2)==8);

figure;
plot(1:length(idx2plot), distance(idx2plot),'bs'); grid on
ha = gca;

set(ha,'XTick',1:length(idx2plot));
set(ha,'XTickLabel',leg4plot(idx2plot));
xlim([0.5 length(idx2plot)+0.5])

ylabel(sprintf('Euclidean distance between \ntreatments (%.0fD-space) ',dim2account))
xlabel('Pair labels')

dist2analyse = distance(idx2plot);
% SRT05 = dist2analyse + [0.2059; 0.1737; 0.0793; 0.2376; 0.0086; 0.1097; 0.0954]; % 0.25*rand(size(dist2analyse));
SRT05 = dist2analyse + 0.35*rand(size(dist2analyse));
SRT09 = dist2analyse + [0.1094; 0.1915; 0.1930; 0.0315; 0.1941; 0.1914; 0.0971]; % 0.2*rand(size(dist2analyse));
SRT10 = dist2analyse + [0.8003; 0.1419; 0.4218; 0.9157; 0.7922; 0.9595; 0.6557]; %1.0*rand(size(dist2analyse));

[r05 p05 l u] = corrcoef(dist2analyse,SRT05);
[r05p pp] = corr(dist2analyse,SRT05,'type','pearson');
rp05 = [r05(1,2) p05(1,2) l(1,2) u(1,2)];

[r09 p09 l u] = corrcoef(dist2analyse,SRT09);
rp09 = [r09(1,2) p09(1,2) l(1,2) u(1,2)];

[r10 p10 l u] = corrcoef(dist2analyse,SRT10);
rp10 = [r10(1,2) p10(1,2) l(1,2) u(1,2)];

xoffset = 0.005;
figure;
plot(dist2analyse-xoffset,SRT05,'x'), hold on, grid on
plot(dist2analyse+xoffset,SRT09,'rs')

dist2ana = round(100*sort(dist2analyse))/100;
set(gca,'XTick',sort(dist2analyse))
set(gca,'XTickLabel',dist2ana)
xlim(minmax(dist2analyse')+[-0.05 +0.05])

[rp05; rp09; rp10]

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

function d = il_euclidean_dist(x,y)

d = sqrt( sum((x-y).^2) );

