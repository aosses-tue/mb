function r20160512_van_Dorp
% function r20160512_van_Dorp
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/05/2016
% Last update on: 12/05/2016 
% Last use on   : 12/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('D:\Databases\dir04-Psychoacoustics\RAA-Jasper-van-Dorp\01-Software\matlab\');
dir = 'D:\Documenten-TUe\06-Lectures-TUe\HIT-2015-2016\Q4\Sounds\';
                                            % sClar     sRev
file = {'T3-sample-Klavier.wav', ...        % 0.5524    0.1504
        'T3-sample-R5-Rij24-Stoel-6.wav'};  % 0.4527    0.1095

for i = 1:length(file)
    [par psi] = raa([dir file{i}]);
    sClar(i,1) = par.sClar;
    sRev(i,1) = par.sRev;  
    
    PsiL = psi.PsiL;
    PsiR = psi.PsiR;
    
    %%%
    Psimin = 7.49e-3; % to split auditory stream for a time > Tmin
    Psimin_dip = -1.33e-3;
    Tmin = 63.1e-3; % 63.1 [ms]

    % for clarity
    a = 2.76;
    b = 2.84;
    
    fs = 16000; % prior knowledge
    [PsiL_dir PsiL_rev] = il_split_stream(PsiL,fs,Psimin,Psimin_dip,Tmin);
    [PsiR_dir PsiR_rev] = il_split_stream(PsiR,fs,Psimin,Psimin_dip,Tmin);

    PsiDIR = sqrt(PsiL_dir.^2 + PsiR_dir.^2);
    PsiREV = sqrt(PsiL_rev.^2 + PsiR_rev.^2);
    K = psi.numBands;
    N = psi.numSamples;
    LDIR = 1/(N*K)*sum(sum(PsiDIR));
    LREV = 1/(N*K)*sum(sum(PsiREV));

    PCLA1(i,1) = LDIR/LREV;
    sPCLA1(i,1) = 1/(1+exp(-a*(PCLA1(i,1)-b)));
    
    Ni  = 4000;
    Nf  = size(PsiL,1);
    N   = Nf-Ni + 1;
    PsiDIR = PsiDIR(Ni:Nf,:);
    PsiREV = PsiREV(Ni:Nf,:);
    LDIR = 1/(N*K)*sum(sum(PsiDIR));
    LREV = 1/(N*K)*sum(sum(PsiREV));

    PCLA2(i,1) = LDIR/LREV;
    sPCLA2(i,1) = 1/(1+exp(-a*(PCLA2(i,1)-b)));

    Ni  = 4000;
    Nf  = 4000+round(5*fs)-1;
    N   = Nf-Ni + 1;
    PsiDIR = PsiDIR(Ni:Nf,:);
    PsiREV = PsiREV(Ni:Nf,:);
    LDIR = 1/(N*K)*sum(sum(PsiDIR));
    LREV = 1/(N*K)*sum(sum(PsiREV));

    PCLA3(i,1) = LDIR/LREV;
    sPCLA3(i,1) = 1/(1+exp(-a*(PCLA3(i,1)-b)));
    
end

disp('')
    % How to calculate sClar, sRev from psi values:
    %   - PCLA = LDIR/LREV; LDIR: Eq 14; LREV: Eq 11
    %   - PASW: Eq 18
    %   - PLEV: Eq 23
    %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outDir outRev] = il_split_stream(Stream,fs,mumin,mumin_dip,Tmin)

outDir  = zeros(size(Stream));
outRev  = zeros(size(Stream));
% iblocks = zeros(size(Stream));

N = size(Stream,1);

for i = 1:size(Stream,2)
    
    Lpsi = 1/N*sum(Stream(:,i));
    Psimin     =     mumin*Lpsi;
    Psimin_dip = mumin_dip*Lpsi;
    
    idx = find(Stream(:,i) >= Psimin);
    idxDir_above = il_detect_segment(idx,fs,Tmin);
    idx = find(Stream(:,i) < Psimin_dip);
    idxDir_below = il_detect_segment(idx,fs,Tmin);
    idxDir = sort([idxDir_above idxDir_below]);
        
    idxRev = 1:size(Stream(:,i),1);
    idxRev(idxDir) = [];
    
    outDir(idxDir,i) = Stream(idxDir,i);
    outRev(idxRev,i) = Stream(idxRev,i);
    
%     figure;
%     plot(outDir(:,i)); hold on
%     plot(outRev(:,i),'r');
    
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idxDir = il_detect_segment(idx,fs,Tmin)

% iblocks = zeros( 1,size(Stream(:,i),1) );
iblocks(idx) = il_detect_blocks(idx);
nBlocks = max(iblocks);
for k = 1:nBlocks
    idxL(k) = find(iblocks==k,1,'first');
    idxU(k) = find(iblocks==k,1,'last');
    if k == 321
        disp('')
    end
end
dur = (idxU-idxL)/fs;
iBlocks2use = find(dur(:) >= Tmin);

idxDir = [];
for k = 1:length(iBlocks2use)
    idx2usetmp = find(iblocks==iBlocks2use(k));
    idxDir = [idxDir idx2usetmp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iblocks = il_detect_blocks(idx);

nBlock = 1;
nConsecutive  =diff(idx);

for i = 1:length(idx)-1
    if nConsecutive(i) == 1
        iblocks(i) = nBlock;
    else
        nBlock = nBlock + 1;
        iblocks(i) = nBlock;
    end
end
iblocks(end+1) = iblocks(end);

