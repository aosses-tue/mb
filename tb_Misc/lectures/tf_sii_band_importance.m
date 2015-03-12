function tf_sii_band_importance
% function tf_sii_band_importance
%
% 1. Description:
%       Band Importance as a function of one-third OB
% 
%       Tested: yes
% 
% Created by Tom Francart
% Edited by Alejandro Osses
% Last use: 10/03/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Copy-paste from scripts Rhebergen:


% Critical band
% cf=[150 250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 ...
%     2500 2900 3400 4000 4800 5800 7000 8500];
% default=        [0.0103, 0.0261, 0.0419, 0.0577, 0.0577, 0.0577, 0.0577,0.0577, 0.0577,0.0577, 0.0577,0.0577, 0.0577,0.0577, 0.0577,0.0577, 0.0577, 0.0460, 0.0343, 0.0226, 0.0110];
% SPIN=           [0.0130, 0.0478, 0.0451, 0.0470, 0.0523, 0.0591, 0.0591, 0.0503, 0.0503, 0.0556, 0.0699, 0.0625, 0.0602, 0.0684, 0.0638, 0.0605, 0.0534, 0.0394, 0.0291, 0.0132, 0.0];

% 1/3 octave band 
f = [160 200 250 315 400 500 630 800 1000 1250 1600 2000, 2500 3150 4000 5000 6300 8000];
 
% columns of BIArr (Band Importance array):
%       1:	Average speech as specified in Table 3
%		2:	various nonsense syllable tests where most English phonemes occur equally often
%		3:	CID-22
%		4:	NU6
%		5:	Diagnostic Rhyme test
%		6:	short passages of easy reading material
%		7:	SPIN

BIArr= [0.0083	0		0.0365	0.0168	0		0.0114	0
		0.0095	0		0.0279	0.013	0.024	0.0153	0.0255
		0.015	0.0153	0.0405	0.0211	0.033	0.0179	0.0256
		0.0289	0.0284	0.05	0.0344	0.039	0.0558	0.036
		0.044	0.0363	0.053	0.0517	0.0571	0.0898	0.0362
		0.0578	0.0422	0.0518	0.0737	0.0691	0.0944	0.0514
		0.0653	0.0509	0.0514	0.0658	0.0781	0.0709	0.0616
		0.0711	0.0584	0.0575	0.0644	0.0751	0.066	0.077
		0.0818	0.0667	0.0717	0.0664	0.0781	0.0628	0.0718
		0.0844	0.0774	0.0873	0.0802	0.0811	0.0672	0.0718
		0.0882	0.0893	0.0902	0.0987	0.0961	0.0747	0.1075
		0.0898	0.1104	0.0938	0.1171	0.0901	0.0755	0.0921
		0.0868	0.112	0.0928	0.0932	0.0781	0.082	0.1026
		0.0844	0.0981	0.0678	0.0783	0.0691	0.0808	0.0922
		0.0771	0.0867	0.0498	0.0562	0.048	0.0483	0.0719
		0.0527	0.0728	0.0312	0.0337	0.033	0.0453	0.0461
		0.0364	0.0551	0.0215	0.0177	0.027	0.0274	0.0306
		0.0185	0		0.0253	0.0176	0.024	0.0145	0];



figure
plot( f, BIArr(:,[1 7]), 'LineWidth', 2);
legend('Average speech', 'SPIN');
grid on
xlabel('1/3 octave band centre frequency (Hz)');
ylabel('Band importance');
set(gca, 'XScale', 'log');
ticks=f(1:2:end);
set(gca, 'XTick', ticks);
set(gca, 'XTickLabel', ticks);
ylim([0 0.105])

end