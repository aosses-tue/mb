function definput = arg_modfilterbank(definput)
% function definput = arg_modfilterbank(definput)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 09/10/2015
% Last update on: 09/10/2015 
% Last use on   : 09/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.modfiltertype    = {'modfilterbank','lowpass'};
definput.flags.modfilterlimit   = {'modfilter_onequarter','modfilter_1kHz_limit'};
definput.flags.modfilter_150Hz_LP = {'LP_150_Hz','no_LP_150_Hz'};
definput.flags.att_factor       = {'att_factor','no_att_factor'}; % Attenuation factor applied to mod filters above 10 Hz
definput.flags.resample_intrep  = {'noresample_intrep','resample_intrep'};

definput.groups.dau1997 = { 'modfilterbank', ...
                            'modfilter_1kHz_limit',...
                            'no_LP_150_Hz', ...
                            'no_att_factor'};

definput.groups.dau1997wLP = { 'modfilterbank', ... % not used in any publication (to my knowledge)
                            'modfilter_1kHz_limit',...
                            'LP_150_Hz', ...
                            'no_att_factor'};
                        
definput.groups.derleth2000 = { 'modfilterbank', ...
                            'modfilter_onequarter',...
                            'no_LP_150_Hz', ...
                            'no_att_factor'};


definput.groups.jepsen2008 = {'modfilterbank', ...
                            'modfilter_onequarter',...
                            'LP_150_Hz', ...
                            'att_factor'};
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
