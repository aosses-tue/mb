clear all, clc, close all

List_of_files = { 'beamformerFilters.mat' };  % Required for Recreating_filt
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

Import_txt_from_ARTA('/home/alejandro/Documenten/Meas/20130218_Meas_SP15_FR_Reference/frs_front.txt');
data_front = data;

Import_txt_from_ARTA('/home/alejandro/Documenten/Meas/20130218_Meas_SP15_FR_Reference/frs_rear.txt');
data_rear = data;

Import_txt_from_ARTA('/home/alejandro/Documenten/Meas/20130218_Meas_SP15_FR_Reference/frs_reference.txt');
data_ref = data;

xmax = 8000; % 8 kHz 

f = data_rear(:,1);

mag_front   = 10.^(data_front     / 20);
mag_rear    = 10.^(data_rear(:,2) / 20);
mag_ref     = 10.^(data_ref(:,2)  / 20);

semilogx(   f, To_dB(mag_front), ...
            f, To_dB(mag_rear ), ...
            f, To_dB(mag_ref  ) )
xlim([1e2 xmax]), grid on
legend('front', 'rear', 'reference','Location','Best')

Count_at_1kHz = max( find(f < 1000) );

front2ref   = mag_front ./ mag_ref;
rear2ref    = mag_rear  ./ mag_ref;

front2ref_norm  = front2ref / front2ref(Count_at_1kHz);
rear2ref_norm   = rear2ref  / front2ref(Count_at_1kHz);

front_final = 20*log10(front2ref_norm);
rear_final  = 20*log10(rear2ref_norm);

figure,
semilogx(   f, 20*log10( rear2ref ) ), hold on
        
legend('Division')
title('Rear')

figure, 
semilogx(   f, front_final, ...
            f, rear_final )
title('Microphone Response Measurement for SP15 Microphones')
xlim([1e2 xmax]), grid on
ylim([-20 20])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]') % Normalised to 1kHz respect to the Front Mic.')
legend('front', 'rear', 'Location', 'Best')

difference = num2str( front_final(Count_at_1kHz) - rear_final(Count_at_1kHz));
display(['Sensitivity difference (Front - Rear) =', num2str(difference) ])

r = Recreate_Filt_SP15_Model(4096);
H = freqz(r.rmic_taps,1,4096/2); H = 20*log10( abs(H) );

maxRear = max( 20*log10( rear2ref ) );
To_sum_H = maxRear - max( H );
H = H + To_sum_H;

H_interpol = interp1(f, 20*log10( rear2ref ) ,r.f);
H_interpol = transpose( H_interpol );
H_after_compensation = 10*log10( (10.^(H_interpol/10)).*(10.^(H/10)) );

figure(2)
plot(r.f, H         , 'r'), grid on % Fs was 8 kHz
plot(r.f, H_interpol, 'g', ...
     r.f, H_after_compensation, 'k')
legend('Response rear','Compensation filter in SP15 Model','Interpolation','after','Location','Best')

xlim([1e2 xmax])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end