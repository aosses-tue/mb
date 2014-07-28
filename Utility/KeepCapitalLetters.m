function y = KeepCapitalLetters(s)

Counter     = ismember(s,'A':'Z');
y = s(Counter);
 