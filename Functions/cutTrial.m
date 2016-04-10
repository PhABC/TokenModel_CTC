function Y = cutTrial(X,start,stop)
    
N = size(X,1);
T = size(X,2);

Y = zeros(N,start+stop);

try
    Y = X(:,end-stop-start+1:end);
catch
    Y(:,end-T+1:end) = X;
    Y(:,1:end-T)     = NaN;
end