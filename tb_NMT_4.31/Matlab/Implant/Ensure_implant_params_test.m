function result = Ensure_implant_params_test

% Ensure_implant_params_test: Test Ensure_implant_params, Ensure_CIC3_params. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No params: default to CIC3

p0 = Ensure_implant_params;
Tester(p0.implant.IC, 'CIC3');
Tester(p0.implant.rf_freq, 5e6);
Tester(p0.implant.protocol, 'embedded');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicitly choose CIC3:

p3.implant.IC = 'CIC3';
p3 = Ensure_implant_params(p3);
Tester(p3.implant.IC, 'CIC3');
Tester(p3.implant.rf_freq, 5e6);
Tester(p3.implant.protocol, 'embedded');
Tester(p3.implant.MIN_PERIOD_us, 69.4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicitly choose CIC4:

p4.implant.IC = 'CIC4';
p4 = Ensure_implant_params(p4);
Tester(p4.implant.IC, 'CIC4');
Tester(p4.implant.rf_freq, 5e6);
Tester(p4.implant.protocol, 'condensed');
Tester(p4.implant.MIN_PERIOD_us, 31.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;
