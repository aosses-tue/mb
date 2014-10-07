function wavfilename = create_MC_CI_stimuli(ref_tone, conf, bGenerateTones, bLB, LB_Gains)
% function wavfilename = create_MC_CI_stimuli(ref_tone, conf, bGenerateTones, bLB, LB_Gains)
%
% Description:
%   ref_tone is a structure with 2 fields (e.g. ref_tone.note  = 'Gsh';
%       ref_tone.octave = 2;)
%   instr          - instruments to be emulated in a cell array 
%                    (e.g. instr = {'UW', 'VIOLIN','CLARINET'};)
%   bGenerateTones - if 0 the wav-files won't be generated but the wavfilenames will
%
%   Run create_APEX3_CI_all_experiments for an example about how to use this script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Matthias M.
% Adapted/Modified by Alejandro Osses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs          = 44100;
instr       = {'UW'}; 

if nargin < 2
    conf = [];
end

if nargin < 3
    bGenerateTones = 1;
end

if nargin < 4
    bLB = 1;
end

if nargin == 0
    ref_tone.note = 'Gsh';
    ref_tone.octave = 2;
    display('Used F0: Gsh_2')
end

wavfilename = [];

cur_dir_old = cd;

if ~isunix
    dirBase = 'C:\Documents and Settings\r0366612\Desktop\Meas';
else
    dirBase = '/home/alejandro/Documenten/Meas/Meas';
end
cur_dir_new = [dirBase delim 'Music' delim 'MC_Stimuli' delim];
try
    cd(cur_dir_new);
catch
    display(['Inside ' mfilename ' : ' cur_dir_new 'does not exist.']);
    cur_dir_new = input('Put a valid directory for saving the new wav-files: ');
    cd(cur_dir_new);
end

f0      = round( note2freq(ref_tone) );

if ~isfield(conf,'ints');       conf.ints = {'1', '2', '3', '4'};   end % default interval separation
if ~isfield(conf,'Strategy');   conf.Strategy = '';                 end

ints        = conf.ints;
Strategy    = conf.Strategy;

seqs    = {'R','FA','FL','FAFL','FLFA','RFA','FAR','RFL','FLR'};

Num_tones_per_serie = length( seqs );
Num_series          = Num_tones_per_serie * length(ints);

dur     = 0.25; % in seconds
gap     = 0.05; % in seconds
dBFS    = -20; % Originally: -6

cont    = 0; % For counting generated WAV files
m       = 1; % Always just 1 instrument

for i   = 1:length(seqs)

   for j= 1:length(ints)

       if bLB == 0
           multisine = create_MCI_sequence([seqs{i} ints{j}], fs, f0, dur, gap, bLB, LB_Gains);
           wavfilenames{j,i} = [instr{m} '-' seqs{i} ints{j} '_' num2str(round(f0)) '_Hz'];
           multisine    = eq_rms(multisine, dBFS); % Tones should come calibrated
       else
           multisine = create_MCI_sequence([seqs{i} ints{j}], fs, f0, dur, gap, bLB, LB_Gains);
           wavfilenames{j,i} = [instr{m} '_LB' Strategy '_' seqs{i} ints{j} '_' num2str(round(f0)) '_Hz'];
       end
       
       if bGenerateTones
            wavwrite( multisine, fs, wavfilenames{j,i} );
            cont = cont + 1;
       end
   end
end
wavfilename = [wavfilename; wavfilenames];

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('Wav files created into directory:')
display(cd)
display([num2str(cont),' WAV-files',' at Fs=',num2str(fs),' Hz. Duration of ',num2str(dur),' seconds'])
display('Please move the generated WAV-files to the corresponding folder')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd(cur_dir_old)