function h = Plot_single_electrode_pitch_EH(test)
% h = Plot_single_electrode_pitch_EH(test)
%
% Outputs:
%    h is the handle of the generated figure
% Inputs: 
%    test is used just to check if the file is included in the path.
%    Call this function without input parameters.
%
% % Example:
%    h = function thesis_chapter1fig05_EHPitchMech;
%
% Adapted from thesis_chapter1fig05_EHPitchMech
% Programmed by Matthias, modified by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
   return;
end
FontSize = 20;

params.duration = 0.4;

s1 = single_channel_sequence(17,  100,  0, [5 13], params);
s2 = single_channel_sequence(12,  100,  0, [5 13], params);
s3 = single_channel_sequence(14,   50,  0, [5 13], params);
s4 = single_channel_sequence(14,  100,  0, [5 13], params);
s5 = single_channel_sequence(14, 1000, 10, [5 13], params);
s6 = single_channel_sequence(14, 1000, 30, [5 13], params);

close all;

nr=3;
nc=2;
margins=[65 50 50 65];
figPos=[50 50 1024 768];
pc = PlotConfig(nr,nc,false);
pc.xGap = 10;
pc.yGap = 10;
pc.margins = margins;
pc.makeFigure(figPos);
pc.makeAxesCoordinates();
pc.makeAxes();

xLim = [0 1e3];
yLim = [6.5 11.5];
for i=1:nr
   for j=1:nc
      pc.axesConfig(i,j).ax_xLim = xLim;
      pc.axesConfig(i,j).ax_yLim = yLim;
      pc.axesConfig(i,j).ax_box = 'on';
      pc.axesConfig(i,j).ax_xTick = 100:100:900;
      pc.axesConfig(i,j).ax_yTick = 5:13;
   end
end
pc.axesConfig(1,1).ax_yTickLabel = (23 - pc.axesConfig(1,1).ax_yTick);
pc.axesConfig(2,1).ax_yTickLabel = (23 - pc.axesConfig(1,1).ax_yTick);
pc.axesConfig(3,1).ax_yTickLabel = (23 - pc.axesConfig(1,1).ax_yTick);
pc.axesConfig(3,1).ax_xTickLabel = pc.axesConfig(3,1).ax_xTick;
pc.axesConfig(3,2).ax_xTickLabel = pc.axesConfig(3,1).ax_xTick;

ep = ElectrodogramPlotter(pc.axesConfig(1,1), 5, 13, 'k', s1, '');
ep.plot();
ep = ElectrodogramPlotter(pc.axesConfig(1,2), 5, 13, 'k', s2, '');
ep.plot();
ep = ElectrodogramPlotter(pc.axesConfig(2,1), 5, 13, 'k', s3, '');
ep.plot();
ep = ElectrodogramPlotter(pc.axesConfig(2,2), 5, 13, 'k', s4, '');
ep.plot();
ep = ElectrodogramPlotter(pc.axesConfig(3,1), 5, 13, 'k', s5, '');
ep.plot();
ep = ElectrodogramPlotter(pc.axesConfig(3,2), 5, 13, 'k', s6, '');
ep.plot();

yl = YLabelConfig(pc.axesConfig(2,1), 'Electrode', pc.yGap(1), 'left', 'middle', 90, FontSize);
yl.plot();
yl = YLabelConfig(pc.axesConfig(1,2), 'Place', pc.yGap(1), 'right', 'middle', 90, FontSize);
yl.plot();
yl = YLabelConfig(pc.axesConfig(2,2), 'Rate', pc.yGap(1), 'right', 'middle', 90, FontSize);
yl.plot();
yl = YLabelConfig(pc.axesConfig(3,2), 'Envelope', pc.yGap(1), 'right', 'middle', 90, FontSize);
yl.plot();
xl = XLabelConfig(pc.axesConfig(3,1), 'Time (ms)', pc.xGap(1), 'right', 'bottom', 0, FontSize);
xl.plot();

h = gcf;
set(h, 'PaperPositionMode','auto')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end