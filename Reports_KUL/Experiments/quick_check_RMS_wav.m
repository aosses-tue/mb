function difference = quick_check_RMS_wav(filename, refdBFS)
% function quick_check_RMS_wav(filename, refdBFS)
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    refdBFS = -20;
end

if nargin == 0
    directory = uigetdir('*','PR: Select an APEX experiment result file');
    filename1 = {   'UW_LB_ACE_104_Hz.wav'; ...
                    'UW_LB_ACE_131_Hz.wav'; ...
                    'UW_LB_ACE_147_Hz.wav'; ...
                    'UW_LB_ACE_208_Hz.wav'; ...
                    'UW_LB_ACE_294_Hz.wav'; ...
                    'UW_LB_F0m_104_Hz.wav'; ...
                    'UW_LB_F0m_131_Hz.wav'; ...
                    'UW_LB_F0m_147_Hz.wav'; ...
                    'UW_LB_F0m_208_Hz.wav'; ...
                    'UW_LB_F0m_294_Hz.wav'};    % /home/alejandro/Documenten/Meas/Meas/Music/PR_Stimuli/UW_LB_ACE_156_Hz.wav
    %filename1 = [directory, filename_part1];
else
    filename1 = filename;
end

if size(filename1,1) == 1 || ischar(filename1)
    directory = '';
    filename = filename1;
    filename1 = {filename};
end

for i = 1:length(filename1)
    try
        x = wavread(filename1{i});
    catch
        try
            directory = 'C:\Documents and Settings\r0366612\Desktop\Meas\Music\PR_Stimuli_new\';
            x = wavread([directory filename1{i}]);
        catch
            directory = 'C:\Documents and Settings\r0366612\Desktop\Meas\Music\PR_Stimuli\';
            x = wavread([directory filename1{i}]);
        end
    end

    difference = rmsdb(x)-refdBFS;

    if abs(difference) < 0.1
        difference = 0;
    end
    
    
    % disp([mfilename '.m: the following is the directory read - ' directory])

    if nargout == 0
        disp(['Directory: ' directory])
        disp([filename1{i} ': delta dB ref. ' num2str(refdBFS) ' = ' num2str(difference) ' dB']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end