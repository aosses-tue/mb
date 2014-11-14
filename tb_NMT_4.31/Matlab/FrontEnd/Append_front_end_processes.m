function p = Append_front_end_processes(p)
% function p = Append_front_end_processes(p)
%
% Append front end processes to a param struct.
%
% Appends the following processes:
%       SP12 (Freedom SP)               SP15 (CP810)
%       -------------------------------------------------------------------
%       Wav_proc                        Wav_proc
%       Normalize_wav_proc* (ExpORL)    Normalize_wav_proc* (ExpORL)
%       (Freedom_microphone_proc**)     -
%       FreedomMicResp_proc (ExpORL)    -
%       -                               Short_term_averager_proc (ExpORL)
%       Input_scaling_proc*             -
%       -                               CP810MicResp_proc (ExpoORL)
%       -------------------------------------------------------------------
%  (* ) Proceses marked with * are added when the boolean variable related 
%       to them is equal to 1
%  (**) Freedom_microphone_proc was replaced by FreedomMicResp_proc 
%       (programmed by Matthias Milczynski, ExpORL), since the original 
%       process was not released to us officially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%    $Revision: #1 $
%        $Date: 2008/03/04 $
%      Authors: Brett Swanson
%
% Edited by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = Ensure_field(p,'normalize',0);
p = Append_process(p, 'Wav_proc');

if p.normalize == 1
    p = Append_process(p, 'Normalize_wav_proc');
end

if isfield(p, 'source');

	switch p.source
        case {'HS8', 'Headset'}
            p = Append_process(p, 'HS8_microphone_proc');
            %%% Added by Sean Lineaweaver 27.Oct.06
        case {'Freedom','SP12'}
            % p = Append_process(p, 'Freedom_microphone_proc');
            p = Append_process(p, 'FreedomMicResp_proc');
        case {'CP810','SP15'} % 
            p = Append_process(p, 'Short_term_averager_proc');
            %%% Added by Alejandro Osses 07.May.13
	% Other sources could be modelled here, e.g. lapel microphone.
	end	

	if ~isequal(p.source, 'Digital')
		p = Append_process(p, 'SPrint_front_end_proc');
		if p.pre_emphasis
			p = Append_process(p, 'Pre_emphasis_proc');
		end
    end
	
    if ~isequal(p.source,'CP810') && ~isequal(p.source,'SP15') % All the cases except the SP15 BTE
        if p.input_scaling
            p = Append_process(p, 'Input_scaling_proc');
        end
    else
        p = Append_process(p, 'CP810MicResp_proc');
        %%% Added by Alejandro Osses 07.May.13
    end
end

end