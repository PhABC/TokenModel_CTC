%Functions for cisekModel1b

%for interactions within PMd, replace kernel K with a NL transfer:
function y = fct(x)
y = (1 ./ (0.3 + exp(-4 * (x-1.3)))) + 0.3;
end