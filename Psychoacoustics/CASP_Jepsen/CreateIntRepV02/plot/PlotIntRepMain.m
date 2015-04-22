%--------------------------------------------------------------------------
% PlotIntRepMain.m
%--------------------------------------------------------------------------
%   plots the output of the pre-processing stages of pemo or casp of an
%   arbitrary input signal and/or the internal representation
%
% usage
%   PlotIntRepMain(IntRep, BM, MF, fs, mode)
%
% input
%   IntRep          : structure containing internal representation
%   BM              : basilar filterbank coefficients & stuff
%   MF              : modulation filterbank coefficients & stuff
%   fs              : sample rate of the input signal
%   mode            : plot modus, possible options are:
%                       'input'       : input signal
%                       'basilar'     : output basilar-membrane filerbank
%                       'hc'          : output hair cell transformation stage
%                       'adapt'       : output adaptation loops
%                       'IntRep'      : internal representation
%                                         - IntRep (Dau et al., 1997a, Fig.6, p.2898),
%                                           three-dimensional: one auditory filter
%                                           (on-freq), model units across modulation
%                                           filters across time
%                                         - IntRep (Dau et al., 1996a, Fig.7, p.3619),
%                                           two-dimensional: one auditory filter
%                                           (on-freq), one modulation filter, model
%                                           units across time
%                 'IntRepImage'     : energy-distribution of IntRep across auditory
%                                 filters and modulation filters
%                 'all'         : plots all options
%
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function PlotIntRepMain(IntRep, BM, MF, fs, mode)

if nargin < 5, error('not enough input arguments'); end

%% settings
time1 = 0:1/fs:(length(IntRep.BM) - 1)/fs;
time2 = time1(1:IntRep.resampleFac:end);
freq1 = round(BM.CenterFreqs);
freq2 = round(MF.CenterFreq);
set(0,'DefaultAxesFontsize',11,'DefaultTextFontsize',11);

%% plot
switch mode
    case 'input'
        % plot input signal
        PlotInput(time1,IntRep)
        
    case 'basilar'
        % amplitude output signal basilar-membrane filterbank across auditory filters across time
        PlotBasilar(time1,freq1,IntRep)
        
    case 'basilarrms'
        % energy distribution (rms across time) output signal basilar-membrane filterbank across auditory filters
        PlotBasilarRms(freq1,IntRep)
        
    case 'hc'
        % amplitude output signal hair cell transformation stage across
        PlotHc(time1,freq1,IntRep)
        
    case 'adapt'
        % output signal in model units after adaptation loops across auditory filters across timefigure
        PlotAdapt(time1,freq1,IntRep)
        
    case 'IntRep'
        % internal representation (several options)
        PlotIntRep(time2,freq1,freq2,IntRep)
        
    case 'IntRepImage'
        % energy distribution internal representation across auditory filters across modulation filters
        PlotIntRepImage(freq1,freq2,IntRep)
        
    case 'all'
        % in demo modus: plots all stages and internal representation
        PlotInput(time1,IntRep)
        PlotBasilar(time1,freq1,IntRep)
        PlotHc(time1,freq1,IntRep)
        PlotAdapt(time1,freq1,IntRep)
        PlotIntRep(time2,freq1,freq2,IntRep)
        if length(freq1)> 1 && length(freq2) > 1
        PlotIntRepImage(freq1,freq2,IntRep)
        else
        end
    otherwise
        error('create IntRep: illegal plot modus')
end

%eof
