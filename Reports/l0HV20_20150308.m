function l0HV20_20150308
% function l0HV20_20150308
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 08/03/2015
% Last update on: 08/03/2015
% Last use on   : 08/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
dir = [Get_TUe_paths('MATLAB') 'Reports' delim];
% file1 = [dir mfilename delim 'Diff_gear-grayscale.jpg'];
file1 = [dir mfilename delim 'Effects_gray.png'];
file2 = [dir mfilename delim 'Effects_color.png'];

%% 0.0 Learning a little bit
% Color images are stored in an RGB array (3D) while grayscaled images are 
% in a 1D array. 

im1 = imread(file1);
[r2 g2 b2] = imread(file2);

im1 = rgb2gray(im1);

figure;
subplot(1,2,1)
imshow(im1)

subplot(1,2,2)
imshow(r2)

%% Experiment A

M = 51;
x = -(M/2-1):M/2-1;
y = -(M/2-1):M/2-1;

a = 2;
b = 1;

expSigma = [-1:5];

for i = 1:length(expSigma)
    
    sigma_a = 2.^expSigma(i); % 7 conditions
    sigma_b = sqrt(a/b)*sigma_a;

    m = a*exp(-(x.^2+y.^2)/(2*sigma_a^2)) - b*exp(-(x.^2+y.^2)/(2*sigma_b^2));

    exp2eval = sprintf('im2%.0f = imfilter(im1,m,''conv'');',i);
    eval(exp2eval);
    
end

figure;
subplot(2,2,1)
imshow(im1);
title('Original')

subplot(2,2,2)
imshow(im21);
title('1')

subplot(2,2,3)
imshow(im22);
title('2')

subplot(2,2,4)
imshow(im23);
title('3')

figure;
subplot(2,2,1)
imshow(im24);
title('4')

subplot(2,2,2)
imshow(im25);
title('5')

subplot(2,2,3)
imshow(im26);
title('6')

subplot(2,2,4)
imshow(im27);
title('7')

%% Experiment B

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
