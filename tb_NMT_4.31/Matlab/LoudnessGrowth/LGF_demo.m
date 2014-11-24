function LGF_demo(options)

% LGF_demo: Plot Loudness Growth Function, for various parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
	options = '';
end

B = 0:3:15;
Q = 20:10:40;

xmax = 200;
xx = 0:xmax;	% X axis vector.

xlabel_str = 'Filter envelope amplitude';
ylabel_str = 'Output magnitude';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot typical values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
p.sat_level  = 150;
p.base_level = 4;
p.lgf_Q = 20;
p = LGF_proc(p);

figure;
plot(xx, LGF_proc(p, xx));
axis([0 xmax 0 1]);

% Now show meaning of Q:
v = p.sat_level/sqrt(10);	% Input level 10 dB down from saturation point
m = LGF_proc(p, v);
line([v 0; v v], [0 m; 1 m], 'LineStyle', ':');

s = sprintf('Loudness Growth Function: Base Level = %d, Q = %d', p.base_level, p.lgf_Q);
title(s);
Window_title(s);
xlabel(xlabel_str);
ylabel(ylabel_str);

if isequal(options, 'print')
	print -dbitmap 'LGF_B4Q20';
	print -depsc2  'LGF_B4Q20';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot varying Base level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
p.sat_level = 150;
p.Q = 20;

figure; hold on;
for n = 1:length(B);
	p.base_level = B(n);
	p = LGF_proc(p);
	plot(xx, LGF_proc(p, xx));
end
axis([0 xmax 0 1]);

x = 30;
p.base_level = B(1);
p = LGF_proc(p);
text(x, LGF_proc(p, x), sprintf('B = %d\\rightarrow ', B(1)), 'horizontalAlignment','right');

x = 20;
p.base_level = B(end);
p = LGF_proc(p);
text(x, LGF_proc(p, x), sprintf('\\leftarrowB = %d', B(end)));

t = 'Loudness Growth Function: Q = 20, Base Level varies';
title(t);
Window_title(t);
xlabel(xlabel_str);
ylabel(ylabel_str);

if isequal(options, 'print')
	print -dbitmap 'LGF_B';
	print -depsc2  'LGF_B';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot varying Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
p.sat_level = 150;
p.base_level = 4;

figure; hold on;
for n = 1:length(Q);
	p.Q = Q(n);
	p = LGF_proc(p);
	plot(xx, LGF_proc(p, xx));
end
axis([0 xmax 0 1]);

x = 60;
p.Q = Q(end);
p = LGF_proc(p);
text(x, LGF_proc(p, x), sprintf('\\leftarrowQ = %d',   Q(end)));

x = 30;
p.Q = Q(1);
p = LGF_proc(p);
text(x, LGF_proc(p, x), sprintf('Q = %d\\rightarrow ', Q(1)), 'horizontalAlignment','right');

t = 'Loudness Growth Function: Base Level = 4, Q varies';
title(t);
Window_title(t);
xlabel(xlabel_str);
ylabel(ylabel_str);

if isequal(options, 'print')
	print -dbitmap 'LGF_Q';
	print -depsc2  'LGF_Q';
end
