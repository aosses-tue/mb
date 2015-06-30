% function y = scsin( mpeak, x, base, sat, doPlot )
%
% Scales sinusoid so sinusoid peak to stay below base and saturation level of
% the LGF to always ensure full modulation depth. Used in the F0mod
% procedure
%
% INPUTS
%   mpeak   : maximun in signal to be modulated
%   x       : 'old modulator'
%   base    : base level
%   sat     : saturation level
%   doPlot  : plots the input and output sinewave
% OUTPUTS
%   y       : scaled sinusoid

function y = scsin( mpeak, x, base, sat )

offset = 0.0001;

pos  = 1;
for i = 2:length(x)
    if( x(i) ~= x(1) )
        pos = i;
        break;
    end
end

px = x - 0.5; % bring to [-0.5, 0.5];

[mins, minpos] = lmin( px(pos:end), 1 );
[maxs, maxpos] = lmax( px(pos:end), 1 );

[minmax, minmax_pos] = min(maxs); % minimal maximum peak
[maxmin, maxmin_pos] = max(mins); % maximal minimum peak

%fprintf(1, 'Maxmin: %f at %f, Minmax: %f at %f \n', maxmin + 0.5, t(minpos(maxmin_pos) + pos - 1), minmax + 0.5, t(maxpos(minmax_pos) + pos - 1));
sat  = ( sat  + offset ) - 0.5;
base = ( base - offset ) - 0.5;

maxQ = abs(minmax/sat);
minQ = abs(maxmin/base);
%fprintf(1, 'MaxQ %f, MinQ %f \n', maxQ, minQ);

if( maxQ < 1 ) 
    if( maxQ <=  minQ )
        %disp('Maxima correction');
        px = px*1/maxQ;
        %fprintf(1, 'Fac: %.4f\n',1/maxQ);
    end
elseif ( minQ < 1)
    %disp('Minima correction');
    px = px*1/minQ;
    %fprintf(1, 'Fac: %.4f\n',1/minQ);
end

y = px + 0.5;   
