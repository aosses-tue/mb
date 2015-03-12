function tf_aided_thresholds
% function tf_aided_thresholds

freqs = [250    500     1000    2000    4000    8000];

audiograms =    [20     30      40      40      60      80  
                 30     40      40      50      60      90 
                 40     40      40      40      40      40
                 60     70      80      80      90      100];

ig=nan(size(audiograms));

for i=1:size(audiograms,1)
    
    ig(i,:) = GetInsertionGainsNALRP(audiograms(i,:), freqs);
    
end

ig