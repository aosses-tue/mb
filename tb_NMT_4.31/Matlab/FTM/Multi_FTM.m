function u = Multi_FTM(p, m, title_str)

% Multi_FTM: Plot multiple FTM images in a single figure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Cochlear Ltd
%   $Header: //dsp/Nucleus/_Releases/NMT_4.30/Matlab/FTM/Multi_FTM.m#1 $
% $DateTime: 2008/03/04 14:27:13 $
%   $Change: 86418 $
%   Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = Examine_FTM(p, m, title_str);

a = Get_axes(u.num_FTMs, 1);

n = [0: u.num_time_samples(1) - 1];
t = n / u.sample_rate;

for n = 1:u.num_FTMs
    a = Get_axes;
	image( max(0, real(u.data{n})),...
              'XData', t,...
		      'CDataMapping','scaled');
		      
	pos = get(gca, 'Position');
	set(gca, 'Position', pos .* [1,1,1,0.95]);

    text(u.duration/20, u.num_bands * 0.8, u.title_str{n}, 'Color', 'white');
	
	if a.row == 1
		xlabel(u.time_label_string);
	else
		set(gca, 'XTickLabel', []);	
	end
	if a.row == ceil(u.num_FTMs/2)
		ylabel('Channel');
	end
end

a = Get_axes('get');
set(a.axes_vec, 'YDir', 'normal',...
                'CLim', [0, u.overall_max_mag],...
                'TickDir', 'out');
u.a = a;
