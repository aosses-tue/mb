function VoD_synchro
% function VoD_synchro
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example. Results are going to be displayed on screen:
%       VoD_synchro;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/07/2014
% Last update on: 25/07/2014 % Update this date manually
% Last use on   : 28/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ti      = [  0.260 0.859 1.464; ... % mode 2
             0.079 0.472 0.864; ... % mode 3
             NaN NaN NaN; ...
             0.101 0.370 0.638];    % mode 5

tcoil   = [  0.684 1.287 1.893; ... % mode 2
             0.369 0.762 1.154; ... % mode 3
             NaN NaN NaN; ...
             0.285 0.552 0.820];    % mode 5

% tcoil = tcoil*(1 - 0.08); % Time correction... coil passes from tcoil +- 8%
         
misc = Get_VoD_params(0);

Tmodel = repmat(misc.Tmodel,1,3);

Delta = tcoil-ti;

Position = (Delta ./ Tmodel)*100;

Position = Delete_NaN_rows(Position);
Position = Position(:);

[m s] = Get_mean(Position);
fprintf('Coil at Mean = %.4f%%; Std = %.4f of mode period\n',m,s);

Trot = misc.Tmodel;
Omegao = 2*pi./Trot;

Angle_deg = ( repmat(Omegao,1,3).*Delta )/pi*180;
Angle_deg = Delete_NaN_rows(Angle_deg);
tmp = Angle_deg(:);
[m s] = Get_mean(tmp);
fprintf('Coil angle at Mean = %.1f deg; Std = %.1f respect to initial pos.\n',m,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end