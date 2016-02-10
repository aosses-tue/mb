function	Roughness(TopLevstr,Lev1str,Lev2str)
% function	Roughness(TopLevstr,Lev1str,Lev2str)
%
% Inputs: 
%       Lev1str = 'Calc','plotresult','params','AM', 'FM'
%       Lev2str =   'fcfmInp' ,'fcfmCalc' ,'fcfmplot' ,
%                   'fcfmInp2','fcfmCalc2','fcfmplot2',
% 
% Example:
%       Roughness; % then follow the instructions displayed in the GUI
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_ind	= 0;
b_ind	= 0;
Bcol	= [.23 .22 .7];

if nargin<1
    close all;
	Hfig	=	figure('Name','Roughness','MenuBar','none',...
   	 	  		 'NumberTitle' ,'off', ...
       			 'Resize'      ,'off', ...
       			 'Colormap'    ,[]);
    TopLevstr='HfdMenu';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%																		%%%
%%%	MAIN MENU                                                           %%%
%%%																		%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(TopLevstr,'HfdMenu')
	clear all;
    b_ind	= 1;
    s_ind	= 0;
    Lftspc	=	25;
	Breedte	=	510;
	Atxt	=	{''	
            'This program is a MATLAB implementation of a psychoacoustic roughness model that was created by Aures and optimized by Daniel and Weber. It can be used to make graphs of roughness as a function of various parameters, or to measure the roughness of an input WAV-sound.'
            ''
            'CHOOSE AN OPTION:'
            '' };
	Alabels =	{	'Create a graph';'Roughness measurement of an input WAV-file' };
	ShowTxt(Atxt)
	Abt1	=	uicontrol('style','pushbutton','string','Create a graph','position',[Lftspc 80 Breedte 30],'fontweight','demi','callback','Roughness(''graph'')');
	Abt2	=	uicontrol('style','pushbutton','string','Roughness measurement of an input WAV-file','position',[Lftspc 45 Breedte 30],'fontweight','demi','callback','Roughness(''inpfile'')');
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%																		%%%
%%% Roughness measurement on input WAV-file                             %%%
%%%																	    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(TopLevstr,'inpfile')
   if nargin<2
      Atxt	=	{  'Now select your input file. Then select the number of frames that has to be analysed.',...
            		'','Adjust SPL* [dB] to : ','Enter Framesize :    ','Enter No. of frames to analyse :   ',...
         			'','','* if disabled a powerlevel of -20 dBFS equals a SPL of 60 dB'};
		ShowTxt(Atxt,'l')
      AdjL	=	uicontrol('style','edit','string','60','position',[300 226 130 20],'fontweight','demi','callback','','enable','off');
      Ninp	=	uicontrol('style','edit','string','8192','position',[280 205 150 20],'fontweight','demi','callback','');
      NoFrinp	=	uicontrol('style','edit','string','1','position',[280 184 150 20],'fontweight','demi','callback','');
      AdjTog	=	uicontrol('style','checkbox','position',[280 228 16 16],'fontweight','demi','callback','Roughness(''Tbox'')');

      %read input signals
		[filename,pathname]=uigetfile('*.wav','Choose Input WAV-file',100,100);
      [source,Fs,Nb]	= wavread([pathname,filename]);
   	SigSize			= max(size(source))
		if Fs<40000
		   error('Fs must be over 40 kHz \n')
		end
    	CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 161 150 20],'fontweight','demi','callback','');
		save TmpMat.mat
   	set(CNFRM,'callback','Roughness(''inpfile'',''Calc'')');   
   elseif strcmp(Lev1str,'Calc')
      load TmpMat.mat
      N	=	2*(round(0.5*(str2num(get(Ninp,'string')))));
      NoFr	=	str2num(get(NoFrinp,'string'));
		if N>SigSize
         ShowTxt({'Input file is shorter than framesize!','','Choose another input file or a lower framesize'})
    		Banaan	=	uicontrol('style','pushbutton','string','BACK','position',[220 160 110 20],'fontweight','demi','callback','Roughness(''inpfile'')');
      else   
         if get(AdjTog,'value')
            SPLdes	=	str2num(get(AdjL,'string'));
            AmpCorr	=	db2amp(SPLdes-83)/rms(source(:,1));
            source	=	AmpCorr*source;
         end
         if bitand(Fs==44100,N==8192)
      	   load InitFs44100N8192.mat;
      	elseif bitand(Fs==40960,N==8192)
      	   load InitFs40960N8192.mat;
      	elseif bitand(Fs==48000,N==8192)
      	   load InitFs48000N8192.mat;
      	elseif bitand(Fs==48000,N==16384)
      	   load InitFs48000N16384.mat;
      	else
	   	  ShowTxt({'','','Initializing'});pause(0.03)
    	     InitAll;
 			end
         Hweights;
			NoFrMax	=	floor(2*SigSize/N)-1;
			if NoFrMax>1
            NoFr	=	round(min([NoFrMax NoFr]));
        		NoFr	=	max([NoFr 1]);
  	   		txt	=	{'frame(s) to analyse :	','','frame(s) to go :		','','The roughness of the input signal is now being calculated','','','(it might take a couple of minutes)'};
			else
		   	NoFr	=	1;
		   	txt	=	{'only 1 frame to analyse','','','The roughness of the input signal is now being calculated','','','','(it might take a couple of minutes)'};
			end
			ShowTxt(txt);Roughness('Prog',1,0);			%show progressbar
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 

     	 	if NoFrMax>1
     	    	NFtxt	=	uicontrol('style','text','string',num2str(NoFr),'position',[345 276 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	  	    	NFTGtxt	=	uicontrol('style','text','string',num2str(NoFr),'position',[345 236 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
			end
   	   pause(0.1)
         LocStep	=	floor(SigSize/(NoFr+1));
         Locs	= 	max([LocStep-N/2 0]):LocStep:SigSize-1;
			qbuiten		=	1:N;
         ProgMax	=	NoFr;					%p-bar max value
         for BigLoop	=	1:NoFr;
            ProgCur	=	BigLoop;			%p-bar current value
            if NoFrMax>1
            	set(NFTGtxt,'string',num2str(NoFr-BigLoop+1))
            	pause(0.02)   
         	end
            load TmpMat.mat
            InputSig	=	source(Locs(BigLoop)+qbuiten,1);	
         	RoughBody ;
            Output(BigLoop)	=	R;
            SPL(BigLoop)	=	mean(rms(InputSig));
            if SPL(BigLoop)>0
               SPL(BigLoop)	=	amp2db(SPL(BigLoop))+83;		% -20 dBFS	<-->	60 dB SPL
            else
   				SPL(BigLoop)=	-400;
            end
			end
      	ShowTxt({'Average SPL is :		dB','','Average roughness is :		 Asper'});
      	Rtxt1	=	uicontrol('style','text','string',num2str(round(10*mean(SPL))/10),'position',[345 276 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
      	Rtxt2	=	uicontrol('style','text','string',num2str(round(100*mean(Output))/100),'position',[345 236 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
      	Rave	=	mean(Output); SPLave	=	mean(SPL);Time	=	SigSize/Fs;
         save TmpMat.mat NoFr Output SPL SigSize Fs Output N Rave SPLave Time filename Locs
      	if NoFr>1
      	    Mbt3	=	uicontrol('style','pushbutton','string','Plot Timepath','position',[385 80 150 30],'fontweight','demi','callback','Roughness(''inpfile'',''plotresult'')');
      	end
	   	Mbt1	=	uicontrol('style','pushbutton','string','Show parameters','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''inpfile'',''params'')');
         s_ind=1;
      end
	elseif strcmp(Lev1str,'plotresult')
      b_ind		=	1;
      load TmpMat.mat
      figure('name','Value of Roughness through time','NumberTitle' ,'off')
      Tijd	=	SigSize/Fs;
      Dtijd	=	(Locs+(N/2))/SigSize*Tijd;
      q		=	1:NoFr+2;
      AvR(q)	=	mean(Output);
      Rval	=	[mean(Output) Output mean(Output)];
      hold on;
      for q = 1:NoFr-1
         PFr	=	plot([Dtijd(q)-(.5*N/Fs) Dtijd(q)+(.5*N/Fs)],[Output(q) Output(q)],'m+-');
      end
      q		=	1:NoFr;
      Tlijn	=	[0 Dtijd(q) Tijd];
      PNo	=	plot(Tlijn,Rval,'k.-',Tlijn,AvR,'b:',[Dtijd(NoFr)-(.5*N/Fs) Dtijd(NoFr)+(.5*N/Fs)],[Output(NoFr) Output(NoFr)],'m+-');
      axis([0 Tijd 0 1.3*max(Rval)]);
      xlabel('time [s]'); ylabel('Roughness R [Asper]');
      title(filename);legend(PNo,'Contemporary Roughness','Average Roughness','Framesize'); 
   elseif strcmp(Lev1str,'params')
      s_ind=1;
      load TmpMat.mat
      ShowTxt({'sample frequency:			Hz',...
            	'frame size:					',...
               'average SPL:				dB',...
               'average roughness			Asper'},'l')
      Rtxt1	=	uicontrol('style','text','string',num2str(Fs),'position',[320 276 65 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
      Rtxt2	=	uicontrol('style','text','string',num2str(N),'position',[320 256 65 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
      Rtxt3	=	uicontrol('style','text','string',num2str(SPLave),'position',[320 236 65 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
      Rtxt4	=	uicontrol('style','text','string',num2str(Rave),'position',[320 216 65 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
		Mbt2	=	uicontrol('style','pushbutton','string','Back to main menu','position',[25 45 150 30],'fontweight','demi','callback','Roughness(''HfdMenu'')');
      if NoFr>1
          Mbt3	=	uicontrol('style','pushbutton','string','Plot Timepath','position',[385 80 150 30],'fontweight','demi','callback','Roughness(''inpfile'',''plotresult'')');
      end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%																	    %%%
%%%				Roughness as a function of various parameters			%%%
%%%  																	%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(TopLevstr,'graph')
	if nargin<2
        ShowTxt({'Create a graph of the dependence of roughness on various parameters of AM- or FM-signals.','','','Choose option:'},'l')   
        Mbt2	=	uicontrol('style','pushbutton','string','Roughness of AM-signals','position',[205 195 150 30],'fontweight','demi','callback','Roughness(''graph'',''AM'')');
        Mbt3	=	uicontrol('style','pushbutton','string','Roughness of FM-signals','position',[205 160 150 30],'fontweight','demi','callback','Roughness(''graph'',''FM'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%	Roughness of AM-signals
  	elseif   strcmp(Lev1str,'AM') 
        if nargin==2
            ShowTxt({'		ROUGHNESS OF AM-SIGNALS','','Dependence on modulation frequency:','','','','Dependence on modulation depth:','','Dependende on Bandwidth (AM noise) :'},'l')
            Btn1	=	uicontrol('style','pushbutton','string','R  =  f(Fc,Fmod)  (1)  ','position',[380 242 130 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmInp'')','HorizontalAlignment','left');
            Btn2	=	uicontrol('style','pushbutton','string','R  =  f(Fc,Fmod)  (2)  ','position',[380 205 130 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmInp2'')','HorizontalAlignment','left');
            Btn3	=	uicontrol('style','pushbutton','string','R  =  f(mod. depth)','position',[380 165 130 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''mdInp'')','HorizontalAlignment','left');
            Btn4	=	uicontrol('style','pushbutton','string','R  =  f(BW)','position',[380 125 130 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''noise'')','HorizontalAlignment','left');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R=f(Fc,Fm) (1) 
        elseif  strcmp(Lev2str,'fcfmInp')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequencies* [Hz]:','Modulation Frequencies* [Hz]:','Modulation Depth:','SPL [dB]:','','* = multiple values possible'},'l')
            FCb	=	uicontrol('style','edit','string','250 500 1000 5000','position',[280 246 250 20],'fontweight','demi','callback','');
            FMb	=	uicontrol('style','edit','string','15 22 33','position',[280 225 250 20],'fontweight','demi','callback','');
            MDb	=	uicontrol('style','edit','string','1','position',[280 204 250 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','70','position',[280 183 250 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb FMb MDb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 150 150 20],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmCalc'')');
        elseif  strcmp(Lev2str,'fcfmCalc')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculation R=f(Fc,Fm) 
            FigDat	=	get(gcf,'UserData');
            CFs		=	sort(str2num(get(FigDat(1),'string')));
            MFs		=	sort(str2num(get(FigDat(2),'string')));
            MDpth		=	str2num(get(FigDat(3),'string'));
            SLev		=	str2num(get(FigDat(4),'string'));
            lenCF		=	length(CFs);
            lenMF		=	length(MFs);
            Rmat		=	zeros(lenCF,lenMF);
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
            Roughness('Prog',1,0);ProgMax	=	lenCF*lenMF; %p-bar
            for cntno1	=	1:lenCF
                for	cntno2	=	1:lenMF
                    load TmpMat.mat
                    ProgCur	=	lenCF*(cntno1-1) + cntno2;	%p-bar val
                    InputSig	=	createAM(CFs(cntno1),MFs(cntno2),MDpth,SLev);
                    RoughBody;
                    Rmat(cntno1,cntno2)	=	R;
                end               
            end
            ShowTxt({'','','','FINISHED CALCULATING'})
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmplot'')');
            save TmpMat.mat Rmat CFs MFs MDpth SLev lenCF lenMF;
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(Fc,Fm) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'fcfmplot')
            load TmpMat.mat
            figure('name','Roughness Plot','NumberTitle' ,'off','MenuBar','none')
            q	=	1:lenMF;
            for k=1:lenCF 
                subplot(1,lenCF,k);plot(MFs(q),Rmat(k,q),'k.-');grid on;
                axis([min(MFs) max(MFs) 0 1.13*max(max(Rmat))]);
                xlabel('Mod. Freq. [Hz]'); ylabel('Roughness R [Asper]');
                title({'Center Freq. [Hz] =',num2str(CFs(k))}); 
            end
            b_ind	=	1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R=f(Fc,Fm) (2) 
        elseif  strcmp(Lev2str,'fcfmInp2')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Mod. Freq. range [Hz]:		   to','Modulation Depth :','No. of steps :','SPL [dB]:'},'l')
            FCb	=	uicontrol('style','edit','string','1000','position',[280 247 150 20],'fontweight','demi','callback','');
            FMminb	=	uicontrol('style','edit','string','10','position',[280 226 40 20],'fontweight','demi','callback','');
            FMmaxb	=	uicontrol('style','edit','string','160','position',[390 226 40 20],'fontweight','demi','callback','');
            MDb	=	uicontrol('style','edit','string','1','position',[280 205 150 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','5','position',[280 184 150 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','60','position',[280 163 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb FMminb FMmaxb MDb RESb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 130 150 20],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmCalc2'')');
            %%%%%%%%%%%%%%%%% calculation R=f(Fc,Fm) (2) %%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'fcfmCalc2')
            FigDat	=	get(gcf,'UserData');
            CF			=	str2num(get(FigDat(1),'string'));
            MFmin			=	str2num(get(FigDat(2),'string'));
            MFmax		=	str2num(get(FigDat(3),'string'));
            MDepth	=	str2num(get(FigDat(4),'string'));
            Nosteps	=	round(str2num(get(FigDat(5),'string')));
            MFstp		=	MFmin:(MFmax-MFmin)/((Nosteps)-1):MFmax;
            SLev		=	str2num(get(FigDat(6),'string'));
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
			Rarr		=	zeros(1,Nosteps);
            Roughness('Prog',1,0); ProgMax	=	Nosteps;		%p-bar val
            for MFcntr	=	1:Nosteps
                load TmpMat.mat
                InputSig	=	createAM(CF,MFstp(MFcntr),MDepth,SLev);
                ProgCur		=	MFcntr;								%p-bar val
                RoughBody;
                Rarr(MFcntr)	=	R;
            end
            ShowTxt({'','','','FINISHED CALCULATING'})
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''fcfmplot2'')');
            save TmpMat.mat CF MDepth MFstp SLev Rarr Nosteps MFmin MFmax;
            s_ind	=	1;         
            %%%%%%%%%%%%%%%%% plot R=f(Fc,Fm) (2) %%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'fcfmplot2')
            load TmpMat.mat
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(MFstp,Rarr,'k.-');grid on;
            axis([MFmin MFmax 0 1.13*max(Rarr)]);
            xlabel('Modulation Frequency [Hz]'); ylabel('Roughness R [Asper]');
            title('Roughness as a function of modulation frequency'); 
            b_ind	=	1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R=f(mod. depth) 
        elseif  strcmp(Lev2str,'mdInp')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Modulation Frequency [Hz]:','Modulation Depth Range:		   to','No. of steps :','SPL [dB]:'},'l')
            FCb	=	uicontrol('style','edit','string','1000','position',[280 247 150 20],'fontweight','demi','callback','');
            FMb	=	uicontrol('style','edit','string','70','position',[280 226 150 20],'fontweight','demi','callback','');
            MDminb	=	uicontrol('style','edit','string','0','position',[280 205 40 20],'fontweight','demi','callback','');
            MDmaxb	=	uicontrol('style','edit','string','1.2','position',[390 205 40 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','5','position',[280 184 150 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','70','position',[280 163 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb FMb MDminb MDmaxb RESb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 130 150 20],'fontweight','demi','callback','Roughness(''graph'',''AM'',''mdCalc'')');
            %%%%%%%%%%%%%%%%% calculation R=f(mod. depth) %%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'mdCalc')
            FigDat	=	get(gcf,'UserData');
            CF			=	str2num(get(FigDat(1),'string'));
            MF			=	str2num(get(FigDat(2),'string'));
            MDmin		=	str2num(get(FigDat(3),'string'));
            MDmax		=	str2num(get(FigDat(4),'string'));
            Nosteps	=	round(str2num(get(FigDat(5),'string')));
            MDstp		=	MDmin:(MDmax-MDmin)/((Nosteps)-1):MDmax;
            SLev		=	str2num(get(FigDat(6),'string'));
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
			Rarr		=	zeros(1,Nosteps);
            Roughness('Prog',1,0); ProgMax	=	Nosteps;		%p-bar val
            for MDcntr	=	1:Nosteps
                load TmpMat.mat
                InputSig	=	createAM(CF,MF,MDstp(MDcntr),SLev);
                ProgCur		=	MDcntr;								%p-bar val
                RoughBody;
                Rarr(MDcntr)	=	R;
            end
            ShowTxt({'','','','FINISHED CALCULATING'})
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''mdplot'')');
            save TmpMat.mat CF MF MDstp SLev Rarr Nosteps MDmin MDmax;
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(mod. depth) %%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'mdplot')
            load TmpMat.mat
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(MDstp,Rarr,'k.-');grid on;
            axis([MDmin MDmax 0 1.13*max(Rarr)]);
            xlabel('Modulation depth'); ylabel('Roughness R [Asper]');
            title('Roughness as a function of modulation depth'); 
            b_ind	=	1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R=f(BWnoise) 
        elseif  strcmp(Lev2str,'noise')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Bandwidth range [Hz]:		      to','Modulation Frequency :','Modulation Depth :','No. of steps :','SPL [dB]:'},'l')
            FCb     =	uicontrol('style','edit','string','1000','position',[280 250 150 20],'fontweight','demi','callback','');
            BWminb	=	uicontrol('style','edit','string','10','position',[280 229 60 20],'fontweight','demi','callback','');
            BWmaxb	=	uicontrol('style','edit','string','2000','position',[370 229 60 20],'fontweight','demi','callback','');
            MFb     =	uicontrol('style','edit','string','70','position',[280 208 150 20],'fontweight','demi','callback','');
            MDb     =	uicontrol('style','edit','string','0','position',[280 187 150 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','7','position',[280 166 150 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','70','position',[280 145 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb BWminb BWmaxb MFb MDb RESb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 122 150 20],'fontweight','demi','callback','Roughness(''graph'',''AM'',''NoiseCalc'')');
        %%%%%%%%%%%%%%%%% calculation R=f(BWnoise) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
		elseif  strcmp(Lev2str,'NoiseCalc')
            FigDat	=	get(gcf,'UserData');
            CF		=	str2num(get(FigDat(1),'string'));
            BWmin	=	str2num(get(FigDat(2),'string'));
            BWmax	=	str2num(get(FigDat(3),'string'));
            MF		=	str2num(get(FigDat(4),'string'));
            MDepth	=	str2num(get(FigDat(5),'string'));
            Nosteps	=	round(str2num(get(FigDat(6),'string')));
            StpSiz	=  log10(BWmax/BWmin)/(Nosteps-1);
            for k=1:Nosteps
                BWstp(k)		=	BWmin*10^((k-1)*StpSiz);
            end
            SLev		=	str2num(get(FigDat(7),'string'));
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
			Rarr		=	zeros(1,Nosteps);
            Roughness('Prog',1,0); ProgMax	=	Nosteps; % progress-bar val
            for BWcntr	=	1:Nosteps
                load TmpMat.mat
                InputSig	=	AMnoise(CF,BWstp(BWcntr),MF,MDepth,SLev);
                ProgCur		=	BWcntr;	% progress-bar value
                RoughBody;
                Rarr(BWcntr)=	R;
            end
            ShowTxt({'','','','FINISHED CALCULATING'})
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''AM'',''noiseplot'')');
            save TmpMat.mat CF MF MDepth SLev Rarr Nosteps BWmin BWmax BWstp;
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(mod. depth) %%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'noiseplot')
            load TmpMat.mat
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(BWstp,Rarr,'k.-');grid on;set(gca,'XScale','log')
            axis([BWmin BWmax 0 1.13*max(Rarr)]);
            xlabel('Bandwidth of noise'); ylabel('Roughness R [Asper]');
            title('Roughness of AM noise as a function of Bandwidth'); 
            b_ind	=	1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%	Roughness of FM-signals                                     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif   strcmp(Lev1str,'FM') 
        if nargin==2
            ShowTxt({'		ROUGHNESS OF FM-SIGNALS','','','Dependence on frequency deviation :','','Dependence on modulation frequency :','','Dependence on SPL :'},'l')
            Btn1	=	uicontrol('style','pushbutton','string','R  =  f(dF)    ','position',[380 222 130 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''BWInp'')');
            Btn2	=	uicontrol('style','pushbutton','string','R  =  f(Fmod)','position',[380 182 130 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''FmodInp'')');
            Btn3	=	uicontrol('style','pushbutton','string','R  =  f(SPL)  ','position',[380 142 130 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''SPLInp'')');
            
            % R=f(dF) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'BWInp')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Modulation Frequency  [Hz]:','Freq. Deviation Range [Hz]:	     to','No. of steps :','SPL [dB]:'},'l')
            FCb     =	uicontrol('style','edit','string','1600','position',[280 247 150 20],'fontweight','demi','callback','');
            MFb     =	uicontrol('style','edit','string','70','position',[280 226 150 20],'fontweight','demi','callback','');
            FDEVminb=	uicontrol('style','edit','string','100','position',[280 205 40 20],'fontweight','demi','callback','');
            FDEVmaxb=	uicontrol('style','edit','string','1000','position',[380 205 50 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','5','position',[280 184 150 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','60','position',[280 163 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb MFb FDEVminb FDEVmaxb RESb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 130 150 20],'fontweight','demi','callback','Roughness(''graph'',''FM'',''BWCalc'')');
            %%%%%%%%%%%%%%%%% calculation R=f(dF) %%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'BWCalc')
            FigDat	=	get(gcf,'UserData');
			CF		=	str2num(get(FigDat(1),'string'));
			MF		=	str2num(get(FigDat(2),'string'));
			dFmin	=	str2num(get(FigDat(3),'string'));
			dFmax	=	str2num(get(FigDat(4),'string'));
			Nosteps	=	str2num(get(FigDat(5),'string'));
			SLev	=	str2num(get(FigDat(6),'string'));
            StpSiz	=  log10(dFmax/dFmin)/(Nosteps-1);
            for k=1:Nosteps
                dFArr(k)		=	dFmin*10^((k-1)*StpSiz);
            end
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
            load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
            Rarr		=	zeros(1,Nosteps);
         
            % Reference signal R0
            Roughness('Prog',1,0);ProgMax	=	Nosteps+1;ProgCur	=1;
            InputSig	=	createAM(CF,MF,1,SLev);
            RoughBody;
            Rzero		=	R/.29;
            % Roughness of FM signals         
            Roughness('Prog',Nosteps+1,1);			%show progress
            for cntr	=	1:Nosteps
                load TmpMat.mat
                ProgCur		=	cntr+1;												%show progress
                InputSig	=	createFM(CF,MF,dFArr(cntr),SLev);
                RoughBody;
                Rarr(cntr)	=	R;
            end
            ShowTxt({'FINISHED CALCULATING','','Reference signal is an AM-tone with:','Cf  		=		Hz','Fmod		= 		Hz','md  		= ','L   		=		dB','R0 = R/0.29 	=		Asper'},'l')
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''dFplot'')');
  	    	CFtxt	=	uicontrol('style','text','string',num2str(CF),'position',[225 217 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	MFtxt	=	uicontrol('style','text','string',num2str(MF),'position',[225 197 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	MDtxt	=	uicontrol('style','text','string','1','position',[225 177 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	Ltxt	=	uicontrol('style','text','string',num2str(SLev),'position',[225 157 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	R0txt	=	uicontrol('style','text','string',num2str(Rzero),'position',[225 137 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
            save TmpMat.mat CF MF dFArr SLev Rarr Nosteps dFmin dFmax Rzero
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(dF) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'dFplot')
            load TmpMat.mat;
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(dFArr,(Rarr/Rzero)*100,'k.-');grid on;
            axis([dFmin dFmax 0 (110*max(Rarr)/Rzero)]);
            xlabel('Frequency deviation [Hz]'); ylabel('R/R0 [%]');set(gca,'XScale','log')
            title('Relative roughness as a function of frequency deviation'); 
            b_ind	=	1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%% 	 R=f(Fmod) 
        elseif  strcmp(Lev2str,'FmodInp')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Frequency Deviation [Hz]:','Modulation Freq. Range [Hz]:	     to','No. of steps :','SPL [dB]:'},'l')
            FCb	=	uicontrol('style','edit','string','1600','position',[280 247 150 20],'fontweight','demi','callback','');
            dFb	=	uicontrol('style','edit','string','800','position',[280 226 150 20],'fontweight','demi','callback','');
            MFminb	=	uicontrol('style','edit','string','1','position',[280 205 40 20],'fontweight','demi','callback','');
            MFmaxb	=	uicontrol('style','edit','string','300','position',[390 205 40 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','5','position',[280 184 150 20],'fontweight','demi','callback','');
            SPLb	=	uicontrol('style','edit','string','60','position',[280 163 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb dFb MFminb MFmaxb RESb SPLb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 130 150 20],'fontweight','demi','callback','Roughness(''graph'',''FM'',''FmodCalc'')');
            %%%%%%%%%%%%%%%%% calculation R=f(Fmod) %%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'FmodCalc')
            FigDat	=	get(gcf,'UserData');
            CF			=	str2num(get(FigDat(1),'string'));
            dF			=	str2num(get(FigDat(2),'string'));
            MFmin		=	max([1 str2num(get(FigDat(3),'string'))]);
            MFmax		=	str2num(get(FigDat(4),'string'));
            Nosteps	=	round(str2num(get(FigDat(5),'string')));
            StpSiz	=  log10(MFmax/MFmin)/(Nosteps-1);
            for k=1:Nosteps
                MFstp(k)		=	MFmin*10^((k-1)*StpSiz);
            end
			SLev		=	str2num(get(FigDat(6),'string'));
            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
            Rarr		=	zeros(1,Nosteps);
            Roughness('Prog',1,0);ProgMax	=	Nosteps+1;ProgCur = 1;
            InputSig	=	createFM(CF,70,dF,SLev);
            RoughBody;
            Rzero		=	R;
            for MDcntr	=	1:Nosteps
                load TmpMat.mat
                ProgCur		=MDcntr+1;												%show progress
                InputSig	=	createFM(CF,MFstp(MDcntr),dF,SLev);
                RoughBody;
                Rarr(MDcntr)	=	100*R/Rzero;
            end
            ShowTxt({'','','','FINISHED CALCULATING'})
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''Fmodplot'')');
            save TmpMat.mat CF dF SLev Rarr Nosteps MFmin MFmax MFstp Rzero
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(Fmod) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'Fmodplot')
            load TmpMat.mat;
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(sort([MFstp 70]),[Rarr(find(MFstp<70)) 100 Rarr(find(MFstp>70))],'k.-');grid on;
			hold on,plot(70,100,'ro');legend('R/R0','R0');
            axis([MFmin MFmax 0 (1.1*max(Rarr))]);
            xlabel('Modulation Frequency [Hz]'); ylabel('R/R0 [%]');set(gca,'XScale','log')
            title('Relative roughness as a function of modulation frequency'); 
            b_ind	=	1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%% 	 R=f(SPL) 			%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'SPLInp')
            ShowTxt({'		ENTER INPUT PARAMETERS','','Center Frequency [Hz]:','Frequency Deviation [Hz]:','Modulation Frequency [Hz]:','SPL Range :			       to','No. of steps :'},'l')
            FCb	=	uicontrol('style','edit','string','1600','position',[280 247 150 20],'fontweight','demi','callback','');
            dFb	=	uicontrol('style','edit','string','800','position',[280 226 150 20],'fontweight','demi','callback','');
            MFb 	=	uicontrol('style','edit','string','70','position',[280 205 150 20],'fontweight','demi','callback','');
            SPLminb	=	uicontrol('style','edit','string','40','position',[280 184 60 20],'fontweight','demi','callback','');
            SPLmaxb	=	uicontrol('style','edit','string','80','position',[370 184 60 20],'fontweight','demi','callback','');
            RESb	=	uicontrol('style','edit','string','5','position',[280 163 150 20],'fontweight','demi','callback','');
            set(gcf,'UserData',[FCb dFb MFb SPLminb SPLmaxb RESb]);   
            CNFRM	=	uicontrol('style','pushbutton','string','CONFIRM','position',[280 130 150 20],'fontweight','demi','callback','Roughness(''graph'',''FM'',''SPLCalc'')');
            %%%%%%%%%%%%%%%%% calculation R=f(SPL) %%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'SPLCalc')
            FigDat	=	get(gcf,'UserData');
            CF			=	str2num(get(FigDat(1),'string'));
            dF			=	str2num(get(FigDat(2),'string'));
            MF			=	max([1 str2num(get(FigDat(3),'string'))]);
            Lmin		=	str2num(get(FigDat(4),'string'));
            Lmax		=	str2num(get(FigDat(5),'string'));
            Nosteps	=	round(str2num(get(FigDat(6),'string')));
			Larr		=	Lmin:(Lmax-Lmin)/(Nosteps-1):Lmax     ;Larr(Nosteps)=Lmax;

            ShowTxt({'','','','CALCULATING ROUGHNESS','','','','(it might take a couple of minutes)'})
            pause(0.01)
			load InitFs40960N8192.mat;
            Hweights;
            save TmpMat.mat Barkno MinBf N2 a0 Bf MinExcdB N50 dFs Cf  N Ntop ei Fei N0 Ntop2 gzi Fs N01 TempIn h0 
            Rarr		=	zeros(1,Nosteps);
         
            % Reference signal R0
            Roughness('Prog',1,0);ProgMax	=	Nosteps+1;ProgCur	=1;
            InputSig	=	createAM(CF,MF,1,60);
            RoughBody;
            Rzero		=	R/.29;
            % Roughness of FM signals         
            Roughness('Prog',Nosteps+1,1);			%show progress
            for cntr	=	1:Nosteps
                load TmpMat.mat
                ProgCur		=	cntr+1;												%show progress
                InputSig	=	createFM(CF,MF,dF,Larr(cntr));
                RoughBody;
                Rarr(cntr)	=	R;
            end
            ShowTxt({'FINISHED CALCULATING','','Reference signal is an AM-tone with:','Cf  		=		Hz','Fmod		= 		Hz','md  		= ','L   		=		dB','R0 = R/0.29 	=		Asper'},'l')
            Mbt3	=	uicontrol('style','pushbutton','string','Plot Result','position',[25 80 150 30],'fontweight','demi','callback','Roughness(''graph'',''FM'',''Lplot'')');
  	    	CFtxt	=	uicontrol('style','text','string',num2str(CF),'position',[225 217 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	MFtxt	=	uicontrol('style','text','string',num2str(MF),'position',[225 197 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	MDtxt	=	uicontrol('style','text','string','1','position',[225 177 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	Ltxt	=	uicontrol('style','text','string','60','position',[225 157 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
  	    	R0txt	=	uicontrol('style','text','string',num2str(Rzero),'position',[225 137 45 30],'fontsize',12,'BackgroundColor',Bcol,'ForegroundColor','white','fontweight','demi','HorizontalAlignment','left');
            save TmpMat.mat CF MF dF Rarr Nosteps Lmin Lmax Rzero Larr
            s_ind	=	1;
            %%%%%%%%%%%%%%%%% plot R=f(SPL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  strcmp(Lev2str,'Lplot')
            load TmpMat.mat;
            figure('name','Roughness Plot','NumberTitle' ,'off')
            plot(Larr,100*Rarr/Rzero,'k.-');grid on;
            axis([Lmin Lmax 0 (110*max(Rarr)/Rzero)]);
            xlabel('L [dB]'); ylabel('R/R0 [%]');
            title('Relative roughness as a function of level L'); 
            b_ind	=	1;
         
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% 	Progress Bar
elseif  strcmp(TopLevstr,'Prog')
    b_ind	=	1;
    Chill	=	get(gcf,'Children');
    Brd	=	min([round((Lev2str/Lev1str)*504) 504]);
    if length(Chill)<15
        Tfr	=	uicontrol('style','frame','position',[25 15 510 12],'BackgroundColor',[0.5 0.5 0.5]);
        if Brd>0
            Tbalk	=	uicontrol('style','frame','position',[29 18 Brd 6],'BackgroundColor',[0.0 0.0 0.5]);
        end
    else
        set(Chill(1),'position',[29 18 Brd 6]);
    end
    pause(0.002);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% 	Save results	  %%%%%%%%%%%%%%%%%%%%%%%%%
elseif  strcmp(TopLevstr,'Save')
    save Tmp2.mat;
    clear all;
    FigDat	=	get(gcf,'UserData');
    load TmpMat.mat
    filname	=	get(FigDat(1),'string');
    clear FigDat;
    save(filname);
    load Tmp2.mat;delete Tmp2.mat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%    	Toggle   	  %%%%%%%%%%%%%%%%%%%%%%%%%
elseif  strcmp(TopLevstr,'Tbox')
    load TmpMat.mat;
    Tval	=	get(AdjTog,'value');
    if Tval
        set(AdjL,'enable','on');
    else
        set(AdjL,'enable','off');
    end      
end

if b_ind<1
    Mbt1	=	uicontrol('style','pushbutton','string','Back to main menu','position',[25 45 150 30],'fontweight','demi','callback','Roughness(''HfdMenu'')');
end
if s_ind==1
    Sfr		=	uicontrol('style','frame','position',[254 16 228 48],'fontweight','demi');
    Stx		=	uicontrol('style','text','string','save data to:','position',[260 20 220 40],'fontweight','demi');
    Snm		=	uicontrol('style','edit','string','output.mat','position',[260 20 150 20],'fontweight','demi','callback','')	 ;
    Sbt		=	uicontrol('style','pushbutton','string','Save','position',[420 20 60 20],'fontweight','demi','callback','Roughness(''Save'');');	 
    set(gcf,'UserData',[Snm]);
end   
