function total = get_num_correct_trials(s,c)

% Pitch Ranking considering 4 comparison tones

total = zeros(4,1);

for i = 1:length(s)
    if c(i) == 1
        switch s{i}(10:11)
            case 'F_'
                total(1) = total(1)+1;
            case 'Fs'
                total(2) = total(2)+1;
            case 'G_'
                total(3) = total(3)+1;
            case 'Gs'
                total(4) = total(4)+1;
        end
	end
    
end