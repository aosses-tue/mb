function temp=a3stimulus_empty(id)
% function a3stimulus_empty(id, datablocks,fixedparameters,variableparameters,simultaneous)
% make simple stimulus with n datablocks in parallel
% datablocks is a cell of datablocks
% parameters is a struct of fixed parameters
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/06/2015
% Last update on: 11/06/2015 % Update this date manually
% Last use on   : 11/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp=['<stimulus id="' id '">' lf];
temp=[temp tb '<datablocks>' lf];
temp=[temp tb '</datablocks>' lf];
temp=[temp '</stimulus>' lf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
