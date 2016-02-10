function [outsig, fc, mfc, extra] = Gammatone_analysis(insig,fs)
% function [outsig, fc, mfc, extra] = Gammatone_analysis(insig,fs)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       dirmain = [Get_TUe_data_paths('lx_Text') 'lx2016-02-04-update-validation-FS' delim];
%       dirout  = [dirmain 'Test-battery' delim];
%       file = [dirout '65-FM_tone-fm_4.00-fc_1500-df_300-SPL_70-w_25-fs_44100-N_88200.wav'];
%       [insig fs] = Wavread(file);
%       [outsig, fc, mfc, extra] = Gammatone_analysis(insig,fs);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 2/6/2014
% Last update: 2/6/2014 % Update this date manually
% Last used: 2/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargin = {};

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank: 
%   insig  = N x  1
%   outsig = N x 31 (if 31 gammatone filters with fc between 80 and 8000 Hz)

if keyvals.fhigh > fs/2
    keyvals.fhigh = fs/2;
end

[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

mfc = [];
extra.outsig_db = rmsdb(outsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
