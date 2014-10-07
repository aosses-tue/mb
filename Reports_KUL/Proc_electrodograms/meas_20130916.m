function meas_20130916
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totaltime = 1600;

xlimACE = 1927;
seq = ConvertSequence_txt2nsb('/home/alejandro/Documenten/Meas/20130913-VU-electrodograms/VU_vrouw_131_SNR_99_ACE_new.out');
vACE = ConvertSequence_nsb2NMT(seq);
subplot(2,1,1)
Plot_sequence_for_GUI(vACE)
xlim([xlimACE-100  xlimACE+totaltime])

seq = ConvertSequence_txt2nsb('/home/alejandro/Documenten/Meas/20130913-VU-electrodograms/VU_vrouw_131_SNR_99_F0mod_new.out');
vF0mod = ConvertSequence_nsb2NMT(seq);
subplot(2,1,2)
Plot_sequence_for_GUI(vF0mod)
xlimF0mod = 2253;
xlim([xlimF0mod-100  xlimF0mod+totaltime])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end