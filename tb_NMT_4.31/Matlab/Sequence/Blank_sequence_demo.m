% Blank_sequence_demo: Demonstration of stimulation blanking process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Tim Neal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a standard ACE map
p = ACE_map;

% Set the blanking period and frequency
p.blank_period  = 0.05;%s
p.blank_freq    = 10;%Hz
p = Append_process(p, 'Blank_sequence_proc');

% Run the strategy code and plot the output
q = Process(p,'asa'); 
Plot_sequence(q, '"ASA" with Blanking');
