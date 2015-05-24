function my_answers
% function my_answers
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Question 3:
h_u = [-10 -8 -5 -2 0];
ih1 = [9.5 7.5 4.4 1.4 0.2];

N = 2;

figure(N);
plot(ih1, h_u,'o-')
grid on
xlabel('Height ih1')
ylabel('Resting potential threhold h_u')

% They are inversely proportional

%% Question 4: Sustained response

h_u2 = [-10 -9 -8 -7   -6 -5   -4  -3   -2 -1  0]; % -6, -4, -3 not stored
ih1pos = [3  2  1  6.2 5.2 4.2 3.2 2.2 1.2 0.4 0.2];

figure(N+1)
plot(ih1pos, h_u2,'o-')
grid on
xlabel('Height ih1 for positive field')
ylabel('Resting potential threhold h_u')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
