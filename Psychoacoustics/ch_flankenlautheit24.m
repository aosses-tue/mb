function [ns,lautheit] = ch_flankenlautheit24(kernlautheit)
% function [ns,lautheit] = ch_flankenlautheit24(kernlautheit)
%
% calculates specific loudness pattern and loudness for all time frames 
% calls flankenlautheit_t24.m to calculate specific loundess for single
% time frames
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000 (new version with comments and examples on 06/01/2007)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 18/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns=zeros(length(kernlautheit(:,1)),240);
for i=1:length(kernlautheit(:,1))
   [ns(i,:), lautheit(i,1)]=ch_flankenlautheit_t24(kernlautheit(i,:));
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end