function r20160404_update_proc_ICRA_APEX
% function r20160404_update_proc_ICRA_APEX
%
% 1. Description:
%       Analysing the intermediate results of the ICRA pilot.
% 
% 2. Stand-alone example:
%       r20160404_update_proc_ICRA_APEX;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 01/04/2016
% Last update on: 01/04/2016 
% Last use on   : 01/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDoPlots = 0;
bDoAnaPedestal = 1;

dir_where = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-5-Multiprocedure-auto\S00-Initial-test-AO' delim];

if bDoPlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = {'piano_multi_result-C2-test-1.apr', 'piano_multi_result-A4-test-1.apr', 'piano_multi_Csh5-test-1.apr'};

h = [];
for i = 1:length(f)

    quick_staircases([dir_where f{i}]);
    h(end+1) = Figure2paperfigureT(gcf,2,2);

    Saveas(h(end),['staircase-' num2str(i)]);

end
end

if bDoAnaPedestal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = [dir_where 'Stimuli-A4' delim 'noise-P1t2_P7t2.wav'];
fp = [dir_where 'Stimuli-A4' delim 'pedestal-P1t2_P7t2-SNR-10-dB.wav'];

[xn fs] = Wavread(fn);
[xp   ] = Wavread(fp);

xp = From_dB(10)*xp(1:length(xn));

[dbn t1] = rmsdb_sec(xn,fs,10e-3);
[dbp t2] = rmsdb_sec(xp,fs,10e-3);

figure;
plot(t1,dbn,t2,dbp);
legend('ICRA','pede')

[iin,ssn] = Get_Loudness_MGB( xn         , fs);
[iii,sss] = Get_Loudness_MGB( xp(1:44100), fs);

figure;
plot(1:length(ssn),ssn,1:length(sss),sss);
legend('ICRA','pede')

disp('')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


