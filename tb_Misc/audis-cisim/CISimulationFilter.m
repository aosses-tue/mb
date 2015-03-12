function B = CISimulationFilter(xelec,fs,b,cochlength);
%Script to calculate the spectrum of resynthesis filters for overlapping
%bands with exponential decay
%
%The filter is a FIR approximation of a desired frequency response.
%The desired frequency response is an exponential decay of amplitude
%along the distance of the cochlea.
%Distance to frequency is determined by Greenwood
%
%   B = CISimulationFilter(xelec,fs,b);
%
%   x_elec = the position of the electrode to be simulated (expressed in mm of insertion from the round window)
%   fs = sampling frequency;

%obtain a number of points along the cochlea where we define the desired amplitude
%We use Npoints points
Npoints = 400;
%the points are equally spaced along the cochlea.
x = [1:(cochlength-2)/Npoints:(cochlength-1)];  %a cochlea of length cochlength with first and last mm dropped
%x = [ [1:xelec-2] [xelec-2:0.05:xelec+2] [xelec+2:34] ]
F = greenwood_x2cf(x)/(fs/2);  %The frequencies corresponding to the points along the cochlea

%An exponential decay of gain along the cochlea
A = exp(-b*abs(x-(cochlength-xelec)));

%only use the frequencies that lie in the possible interval (0 to 1.0)
ind = find(F<1);
if mod(length(ind),2) > 0
    ind = ind(1:end-1);
end

B= fir2(150,[0 F(ind) 1],[0 A(ind) 0]);