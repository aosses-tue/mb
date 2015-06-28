function out = Get_psyparams(param)
% function out = Get_psyparams(param)
%
% 1. Description:
%       Loads fixed variables to be used in the Roughness (R) and Fluctuation
%       Strength (FS) models.
% 
%       -----------------
%               R   FS
%       -----------------
%       a0tab   yes yes
%       Bark    yes yes
%       gr      yes no
%       -----------------
% 
% 2. Stand-alone example:
%       a0tab = Get_psyparams('a0tab');
% 
% 3. Additional info:
%       Tested cross-platform: Yes
% 
%       see also Roughness_offline to use this script
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/11/2014
% Last update on: 12/11/2014 % Update this date manually
% Last use on   : 26/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch param
    case 'a0tab' 
        a0tab =	[ 0	 0
          10	 0
          12	 1.15
          13	 2.31
          14	 3.85
          15	 5.62
          16	 6.92
          16.5	 7.38
          17	 6.92
          18	 4.23
          18.5	 2.31
          19	 0
          20	-1.43
          21	-2.59
          21.5	-3.57
          22	-5.19
        22.5	-7.41
          23	-11.3
          23.5	-20
          24	-40
          25	-130
          26	-999];
      out = a0tab;
      
    case 'Bark'
        Bark = [0     0	   50	 0.5
                1   100	  150	 1.5
                2   200	  250	 2.5
                3   300	  350	 3.5
                4   400	  450	 4.5
                5   510	  570	 5.5
                6   630	  700	 6.5
                7   770	  840	 7.5
                8   920	 1000	 8.5
                9  1080	 1170	 9.5
                10  1270 1370	10.5
                11  1480 1600	11.5
                12  1720 1850	12.5
                13  2000 2150	13.5
                14  2320 2500	14.5
                15  2700 2900	15.5
                16  3150 3400	16.5
                17  3700 4000	17.5
                18  4400 4800	18.5
                19  5300 5800	19.5
                20  6400 7000	20.5
                21  7700 8500	21.5
                22  9500 10500	22.5
                23 12000 13500	23.5
                24 15500 20000	24.5];
        out = Bark;
    case 'gr'
        % BarkNo  0     1   2   3   4   5   6   7   8     9     10
        %	 11     12  13  14  15  16  17  18  19  20  21  22  23  24 
        gr = [  0,1,2.5,4.9,6.5,8,9,10,11,11.5,13,17.5,21,24;
                0,0.35,0.7,0.7,1.1,1.25,1.26,1.18,1.08,1,0.66,0.46,0.38,0.3];
        out = gr;
        
    case 'gr-ERB'
        % BarkNo  0     1   2   3   4   5   6   7   8     9     10
        %	 11     12  13  14  15  16  17  18  19  20  21  22  23  24 
        gr = [  0,1,2.5,4.9,6.5,8,9,10,11,11.5,13,17.5,21,24;
                0,0.35,0.7,0.7,1.1,1.25,1.26,1.18,1.08,1,0.66,0.46,0.38,0.3];
        fb      = audtofreq( gr(1,:) ,'bark');
        erbb    = freqtoaud( fb      ,'erb' );
        gr(1,:) = erbb;
        out = gr;
    case 'HTres' %Make list of minimum excitation (Hearing Treshold)
        HTres= [	0		130
                    0.01    70
                    0.17	60
                    0.8     30
                    1		25
                    1.5     20
                    2		15
                    3.3     10
                    4       8.1
                    5		  6.3
                    6		  5
                    8		  3.5
                    10		  2.5
                    12		  1.7
                    13.3	  0
                    15		 -2.5
                    16		 -4
                    17		 -3.7
                    18		 -1.5
                    19		  1.4
                    20		  3.8
                    21		  5
                    22		  7.5
                    23 	 15
                    24 	 48
                    24.5 	 60
                    25		130];
        out = HTres;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
