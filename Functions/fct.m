%Functions for cisekModel1b

%for interactions within PMd, replace kernel K with a NL transfer:
function y = fct(x,steep)

	 y = ( (1 ./ (0.4 + exp(-steep*(x-55)))) + 0.1) .*10;
end

