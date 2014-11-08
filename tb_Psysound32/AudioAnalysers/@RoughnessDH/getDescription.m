function str = getDescription(obj)
% GETDESCRIPTION Returns a text string of the description

str = {'This is the Daniel and Weber Roughness Algorithm,', ...
    'introducing overlap of 50 percent, according to D. Hermes', ...
    'This version has been changed in June of 2009, to improve ', ...
    'the accuracy of the algorithm. Previously a 60dB tone ', ...
    '100% modulated at 70Hz would output an asper value of 0.91.', ...
    'Clearly this should be a value of 1 asper, and the gain was ', ... 
		'increased to effect this. The paper that documents this model is:'...
		'Daniel, P., & Weber, R. (1997). Psychoacoustical roughness: implementation'...
		'of an optimized model. Acustica(83), 113-123.'};
