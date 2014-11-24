function Replot_sequences(c)

% Replot_sequences: Replot a cell array of sequences.
% Replaces the current sequence plot (if any).
% Intended to be used as a Psychophysics present_stimuli function. 
%
% Replot_sequences(c)
%
% c:   Cell array of pulses sequences.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

top_fig = get(0,'CurrentFigure');

[was_seq_fig, old_seq_fig] = figflag('Sequence', 1);

q = Concatenate_sequences(c);
Save_sequence(q, 'q');
Plot_sequence(q, 'Sequence', [1,22]);

if was_seq_fig
	old_seq_fig = old_seq_fig(end);
	pos = get(old_seq_fig, 'Position');
	close(old_seq_fig);
	set(gcf, 'Position', pos);
end

if ~isempty(top_fig)
	figure(top_fig);
end


