function temp=a3datablock(id,filename,device,channels)
% function result = a3datablock(id,filename,device,channels)
%
% 1. Description:
%       id - string
%       filename - string ('*.wav')
%       device - string ('wavdevice')
%       channels - integer
% 
%       lf=sprintf('\n'); % enter
%       tb=sprintf('\t'); % tab
% 
% 2. Stand-alone example:
%       id = 'rawdata-20-0-0';
%       filename = 'sin2000Hz_sin2000Hz_itd0.wav';
%       device = 'wavdevice';
%       result = a3datablock(id,filename,device);
%       disp(result)
% 
% Comments by Alejandro Osses, ExpORL 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<4)
    channels=0;
end

if (nargin<3)
    device='';
end

lf=sprintf('\n');
tb=sprintf('\t');

temp=['<datablock id="' id '">' lf];
if (~isempty(device))
    temp=[temp tb '<device>' device '</device>' lf];
end
temp=[temp tb '<uri>' xmlfilename(filename) '</uri>' lf];
if (channels)
   temp=[temp wraptag('channels', num2str(channels))]; 
end
temp=[temp '</datablock>' lf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%