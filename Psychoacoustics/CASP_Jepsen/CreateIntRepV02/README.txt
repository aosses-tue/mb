% -------------------------------------------------------------------------
% Calculate and Plot the Output of the Pre-Processing Stages and/or the 
% Internal Representation (IntRep) of your stimulus with PEMO or CASP
% -------------------------------------------------------------------------
% version 1.0
%   23/01/2013, C.T. Iben and S.D. Ewert
% version 2.0
%   29/01/2013, C.T. Iben and S.D. Ewert

%%%%%%%%%%%%
% Content: %
%%%%%%%%%%%%
%
% This small package offers the possibility to calculate and plot the output
% of the pre-processing stages and/or the internal representation  of an 
% arbitrary input signal with the models PEMO (Dau et al., 1997a) and CASP
% (Jepsen et al., 2008). 
% 
% Depending on the ranges of the auditory and modulation filterbank several
% plot options for the internal representations are available, e.g.  it is 
% possible to generate two- and three-dimensional internal representations 
% as those published in Dau et al. (1997a, Fig.6, p.2898) and Dau et al. 
% (1996a, Fig.7, p.3619). Additionally the output stages in the pre-processing
% of both models can be visualized. For further details please browse 
% PlotIntRepMain.m, PlotIntRep.m and PlotIntRepImage.m
% 
%%%%%%%%%
% Demo: %
%%%%%%%%%
%
% For a quick start two example stimuli are attached. Just run DemoMain 
% without any input arguments. In this case the following test stimulus will
% be used:
%
%%%%%%%%%%%%%%%%%%%%
% Test Stimulus 1: %
%%%%%%%%%%%%%%%%%%%%
% 
% A 20 Hz sinusoidal amplitude-modulation with a modulation depth of -6 dB 
% was impressed on a 3 Hz wide running (Gaussian) noise carrier centered at
% 5 kHz; stimulus level: 65 dB; stimulus length: 800 ms with 150 ms raised 
% cosine ramps. (saved as x.mat)
% 
%%%%%%%%%%%%%%%%%%%%
% Test Stimulus 2: %
%%%%%%%%%%%%%%%%%%%%
%
% Logatome duhd from OLLO data base (Wesker et al., 2005), male speaker.
%
% In both cases the following model defaults are set:
% 
%%%%%%%%%%%%%%%%%%%
% Default Values: %
%%%%%%%%%%%%%%%%%%%
%
% default model         : Pemo;
%
% range of the basilar-membrane filterbank:
% fcMin                 :    100 Hz;
% fcMax                 :  10000 Hz;
% ERBden                :     1 ERB;
% 
% range of the modulation filterbank:
% fmMin                 :     0 Hz;
% fmMax                 :  1000 Hz;
%
% default plot modus    : all; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Change Stimulus, Models & Stuff: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To change the stimulus, the model or default values please browse the 
% appropriate files:
% 
% for stimulus, model or plot options  ´    : DemoMain.m  
% for filterbank parameter                  : PemoPreProcCfg.m, CaspPreProcCfg.m
% for filterbank coefficients and stuff     : PemoPreProcInit.m, CaspPreProcInit.m
% for plot details                          : PlotIntRepMain.m, PlotIntRep.m,
%                                             PlotIntRepImage.m
%
% All other model parameters are set as in the published model versions and 
% can be looked up in the corresponding pre-processing files in the casp2008 and 
% pemo1997 folder(i.e. PemoPreProc.m,CaspPreProc.m)
%
%%%%%%%%%%%%%%%%%
% Shared files: %
%%%%%%%%%%%%%%%%%
%
% There is an order <shared> in the package including all the m-files used by
% both models. In order to offer some transparency which specific files are
% required to run the pre-processing of each model properly, empty txt-files 
% sharing the same name with the corresponding m-files are inserted in
% the respective model subfolder. 
%
%
%%%%%%%%%%%%%%%
% References: %
%%%%%%%%%%%%%%%
%
% Dau, T. , Püschel, D. and Kohlrausch, A. (1996): "A quantitative model of the
%     `effective' signal processing in the auditory system: (I). Model structure",
%     J. Acoust. Soc. Am. 99, p. 3615-3622.
% 
% Dau, T., Kollmeier, B. and Kohlrausch, A. (1997):"Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band carriers",
      J. Acoust. Soc. Am. 102, p. 2892-2905.
%
% Jepsen, M. L., Ewert, D. S. and Dau, T. (2008):"A computational model of 
%     human auditory signal processing and perception.", 
%     J. Acoust. Soc. Am. 124, p, 422-438.
%
% Meddis, R., O'Mard, L. P. and Lopez-Poveda E. A. (2001): "A computational 
%     algorithm for computing nonlinear auditory frequency selectivity.",
%     J. Acoust. Soc. Am. 109, p. 2852-2861.      
%
% Püschel, D. (1988): "Prinzipien der zeitlichan Analyse beim Hören," Doctoral Thesis,
%     Universität Göttingen.
% Wesker, T., Meyer, B., Wagener, K., Anemüller, J., Mertins, A. and Kollmeier, B.
% (2005): "Oldenburg Logatome Speech Corpus (OLLO) for Speech Recognition
Experiments with Humans and Machines", in proceedings of Interspeech, p.1273-1276.

%eof

