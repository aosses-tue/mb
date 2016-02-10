function ShowTxt(InpStr,Ali)
% function ShowTxt(InpStr,Ali)
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 13/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
Bcol	=	[.23 .22 .7];
Lftspc	=	25;
Breedte	=	510;

Tpos	=	[Lftspc	355	Breedte	35];
Tfr	=	uicontrol('style','frame','position',Tpos,'BackgroundColor',Bcol);
T		=	uicontrol('style','text','string','CALCULATING PSYCHOACOUSTIC ROUGHNESS','position',Tpos+[4 4 -8 -8],'fontsize',12.5,'fontweight','bold','BackgroundColor',Bcol,'ForegroundColor','white');
if nargin>0
	Apos	=	[Lftspc 123 Breedte 188];         
	Afr	=	uicontrol('style','frame','position',Apos,'BackgroundColor',Bcol);
	A		=	uicontrol('style','text','string',InpStr,'position',Apos+[4 4 -8 -8],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi');
end
if nargin>1
   set(A,'HorizontalAlignment','left');
end