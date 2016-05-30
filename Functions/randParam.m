function Y = randParam(X,std,N)
%% randparam 
% This function will output random values centered around the values 
% provided in X with standart deviation of 'std'. The output will be a 
% matrix Y of dimensions (length(X) by N).

d1 = length(X);

% Row vector condition
if size(X,1)>size(X,2); X = X'; end

% Repmat X for vector operations
X = repmat(X,N,1);

% Random values
R = randn(N,d1)*std;
    
% Gaussian distributed parameters  
Y = X + R;
    
