function Y = binmean(X,nbins)
% Will average inputs of dimensions dim to compress 
% signal in nbins.

dim     = size(X);              
binsize = floor(dim(2)/nbins);  

Y = zeros(dim(1),nbins);

for i = 1:nbins
    range = (i-1)*binsize+1:i*binsize;
    Y(:,i) = nanmean(X(:,range),2);
end
 
