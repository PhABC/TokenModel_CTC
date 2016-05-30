function W = wMat(sd,Ww,sunkR,spars,N,wType)
% Creates a weight matrix based on gaussian distribution.
% The weight matrix is wrapped around

% sd    : Standart deviation of the weight distribution
% Ww    : Weight of the weight matrix
% sunkR : The proportion of sunking for the sunken gaussian weigth distribution
% spars : Random sparsity ratio (1 == all connected)
% N     : number of units
% wType : Type of weigths distribution

%% Distribution types
%	0 : Sunken gaussian distribution
% 	1 : Opposite tuning inhibition
%	2 : Mexicain hat distribution (i.e. near inhibition)
% 	3 : Near and opposite tuning inhibition 
%
% Notes: 
%   +You can change the width of near and opposite inhibition
%	      if you use type 1,2 or 3;
%   +For distribution type 3, you would need to adjust the 
%         strength of both inhibition type in the function 
%         below if you want them to be different.

if nargin ~= 6
	wType  = 0; %Distribution type
end

% For option 1,2,3
IEratio   = 2;    % Strength of inhibition with respect to excitation   
sdNearInh = sd*5; % Standart deviation of near tuning inhibition distribution
sdOppInh  = sd;   % Standart deviation of oppositve tuning inhibition distribution

% Cummulative weight matrix
W = gaussianMat(N,sd);

switch wType     
    case 0
    % Sunking gaussian distribution
        sunk = max(max(W))*sunkR;
        W    = W-sunk; 

    case 1
    % Opposite tuning inhibition    
        oppInhW = circshift(gaussianMat(N,sdOppInh),N/2,2);
        W       = W-circshift(oppInhW,N/2,2);    

        %Normalization & weighting
        W(W>0) = W(W>0)/max(W(W>0))/IEratio;    
        W(W<0) = W(W<0)/abs(min(W(W<0)));    

    case 2
    % Mexicain hat
        nearInhW = gaussianMat(N,sdNearInh);
        W        = W*2-nearInhW;

        %Normalization & weighting
        W(W>0) = W(W>0)/max(W(W>0))/IEratio;    
        W(W<0) = W(W<0)/abs(min(W(W<0)));    

    case 3 
    % Near and opposite tuning inhibition
        oppInhW  = circshift(gaussianMat(N,sdOppInh),N/2,2);
        nearInhW = gaussianMat(N,sdNearInh);
        W        = W*7-oppInhW-3*nearInhW;

        %Normalization & weighting
        W(W>0) = W(W>0)/max(W(W>0))/IEratio;    
        W(W<0) = W(W<0)/abs(min(W(W<0)));    
end

W = W*Ww + randn(size(W))*Ww*0.1; % Adding randomness	
W = W.*100/N;	     % Normalize wrt 100 as parameters were found with 100 N
W(rand(N)>spars) = 0;    % Extra random sparsness as function of tuning curve similarity
% W(W<0) = 0;            % To make race model
    	
function mat = gaussianMat(N,sd)

x = 1:N*2;
all_shifts = linspace(1,N*2,N*2);    
all_g = zeros(N*2,N);
k = 1;

for i = 1:N
    s = all_shifts(k);
    all_g(:,k) = exp((-(x-s).^2)./(2*sd.^2));
    k = k+1;
end

% Excitation connections
mat = all_g(1:N,:)+all_g(N+1:end,:)+... 
       flip(flip(all_g(N+1:end,:),2));
