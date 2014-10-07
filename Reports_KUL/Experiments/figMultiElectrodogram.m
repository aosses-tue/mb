function [h ha] = figMultiElectrodogram(info,s1,s2,s3,s4)
% function [h ha] = figMultiElectrodogram(info,s1,s2,s3,s4)
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
% Original name: thesis_figEHPitchMech
%
% Programmed by Matthias, modified by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ha = [];

yShift = 6; % to fit properly legend in NorthEast postition

info = Ensure_field(info,'xlim',[0 1000]);
info = Ensure_field(info,'ylim',[-1 22+yShift]);
info = Ensure_field(info,'title','');
info = Ensure_field(info,'LocationLegend','NorthEast');
nr=3;
nc=2;
margins=[65 65 65 65];
figPos=[0 0 1024 225*3];
pc = PlotConfig(nr,nc,false);
pc.xGap = 10;
pc.yGap = 10;
pc.margins = margins;
pc.makeFigure(figPos);
pc.makeAxesCoordinates();
pc.makeAxes();

for i=1:nr
   for j=1:nc
      pc.axesConfig(i,j).ax_xLim = info.xlim;
      pc.axesConfig(i,j).ax_box = 'on';
      pc.axesConfig(i,j).ax_xTick = info.xlim(1)+20:20:info.xlim(2)-20;
      if i ~=3
          pc.axesConfig(i,j).ax_yLim = info.ylim;
          pc.axesConfig(i,j).ax_yTick = info.ylim(1)+2:4:info.ylim(2)-yShift;
      else
          pc.axesConfig(i,j).ax_yLim = [50 350];
          pc.axesConfig(i,j).ax_yTick = 75:25:300;
      end
   end
end
pc.axesConfig(1,1).ax_yTickLabel = (23 - pc.axesConfig(1,1).ax_yTick);
pc.axesConfig(2,1).ax_yTickLabel = (23 - pc.axesConfig(1,1).ax_yTick);
pc.axesConfig(3,1).ax_yTickLabel = pc.axesConfig(3,1).ax_yTick;
pc.axesConfig(3,1).ax_xTickLabel = pc.axesConfig(3,1).ax_xTick;
pc.axesConfig(3,2).ax_xTickLabel = pc.axesConfig(3,1).ax_xTick;

ep = ElectrodogramPlotter(pc.axesConfig(1,1), info.ylim(1), info.ylim(2), 'k', s1, '');
ep.plot();
ha(end+1) = gca;
title(info.title)
legend('Sim ACE', 'location', info.LocationLegend);

ep = ElectrodogramPlotter(pc.axesConfig(1,2), info.ylim(1), info.ylim(2), 'k', s2, '');
ep.plot();
ha(end+1) = gca;
legend('NMT ACE', 'location', info.LocationLegend);

ep = ElectrodogramPlotter(pc.axesConfig(2,1), info.ylim(1), info.ylim(2), 'k', s3, '');
ep.plot();
ha(end+1) = gca;
legend('Sim F0mod', 'location', info.LocationLegend);

ep = ElectrodogramPlotter(pc.axesConfig(2,2), info.ylim(1), info.ylim(2), 'k', s4, '');
ep.plot();
ha(end+1) = gca;
legend('NMT F0mod', 'location', info.LocationLegend);

ep = ElectrodogramPlotter(pc.axesConfig(3,1), 75, 350, 'k', s4, '');
ep.plot(); % Non visible plots, just to get handle
ha(end+1) = gca;

ep = ElectrodogramPlotter(pc.axesConfig(3,2), 75, 350, 'k', s4, '');
ep.plot(); % Non visible plots, just to get handle
ha(end+1) = gca;

yl = YLabelConfig(pc.axesConfig(2,1), 'Electrode', pc.yGap(1), 'left', 'middle', 90);
yl.plot();
yl = YLabelConfig(pc.axesConfig(2,1), 'Electrode', pc.yGap(1), 'left', 'middle', 90);
yl.plot();
yl = YLabelConfig(pc.axesConfig(3,1), 'F0 [Hz]', pc.yGap(1), 'left', 'middle', 90);
yl.plot();
xl = XLabelConfig(pc.axesConfig(3,1), 'Time (ms)', pc.xGap(1), 'right', 'bottom', 0);
xl.plot();

h = gcf;
set(gcf, 'PaperPositionMode','auto')

linkaxes(ha(:),'x')
xlim( info.xlim )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end