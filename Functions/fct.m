%Functions for cisekModel1b

%for interactions within PMd, replace kernel K with a NL transfer:
function y = fct(x)
y = ((1 ./ (0.3 + exp(-0.16*(x-55)))) + 0.1).*10;
end

