function auxURL=amtauxdataurl(newURL)
%amtauxdataurl URL of the auxiliary data
%   Usage: url=amtauxdataurl
%          amtauxdataurl(newurl)
%
%   `url=amtauxdataurl` returns the URL of the web address containing
%   auxiliary data.
% 
%   `amtauxdataurl(newurl)` sets the URL of the web address for further calls
%   of `amtauxdataurl`.
%
%   See also: amtauxdatapath amtload

persistent AuxDataURL;

if exist('newURL','var')
  AuxDataURL=newURL;
elseif isempty(AuxDataURL)
    % AuxDataURL=['http://www.sofacoustics.org/data/amt-' amthelp('version') '/auxdata'];
    AuxDataURL=['http://www.sofacoustics.org/data/amt-0.9.7/auxdata'];
end
auxURL=AuxDataURL;
