function C = stackStruct(X)
% Will concatenante same fields in diffenret cells
% Cells in X must contain the same fields and dimensions

ncells = length(X);

C  = X{1};
fn = fieldnames(C); 

for c = 2:ncells
    for f = 1:length(fn)-2
        C.(fn{f}) = vertcat(C.(fn{f}),X{c}.(fn{f}));
    end
end

    