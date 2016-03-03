function Y = cutTrial(X,start,stop)

try
    Y = X(:,end-stop-start:end);
catch
    Y = [];
end