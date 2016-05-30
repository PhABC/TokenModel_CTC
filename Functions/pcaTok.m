function [loads,build,scores,latent] = pcaTok(X,norm,conds,comps) 

sym = 1; 	 % If == 1, it will duplicate cells in a symmetrical fashion
	         % with respect to target
if nargin == 1
    norm  = 0;
    comps = 1:6;
    if sym ; conds = 1:6; else conds = 1:3; end

elseif nargin == 2
    comps = 1:6;
    if sym ; conds = 1:6; else conds = 1:3; end
elseif nargin == 3
    comps = 1:6;
end

nbNets = length(X);

% Stack cells wrt fields
Xnet = cell(nbNets,1);
for net = 1:nbNets;
    Xnet{net} = stackStruct(X{net});  % stackStruct is adapted to this FR format
end

Xfr = stackStruct(Xnet);

Xfr.full   = vertcat(Xfr.PMd,Xfr.M1);
Xfr.nbNets = nbNets; %To feed the variable to function symCells

% Symetrise cells for left direction
if sym
	X2fr = symCells(Xfr);
	X2fr = structfun( @(X) X(:,:,conds), X2fr, ...
                   	'UniformOutput', false);
else
        rmfields = {'Idx','Com','nbTokens','pref','npref','nbNets'};
	Xfr      = rmfield(Xfr,rmfields);

        X2fr = structfun( @(X) X(:,:,conds), Xfr, ...
                   	'UniformOutput', false)

end

% Square root transform
X2fr = structfun( @(X) sqrt(X),X2fr, ... 
                   'UniformOutput', false);

% Normalizing for max FR
if norm
	X2fr = structfun( @(X) X/max(max(max(X),[],3),[],2),X2fr,...
			'UniformOutput', false);
end

fn = fieldnames(X2fr); 
for f = 1:length(fn)
    N = size(X2fr.(fn{f}),1);
    T = size(X2fr.(fn{f}),2);
    C = size(X2fr.(fn{f}),3);

    Xresh.(fn{f}) = reshape(X2fr.(fn{f})(:,5:T-8,:),N,(T-5-7)*C);
    [loads.(fn{f}),scores.(fn{f}),latent.(fn{f})] = pca(Xresh.(fn{f})'); 
end
 
%% Building components
for f = 1:length(fn)
    for comp = comps
        a = 0;
        for cond = 1:length(conds)
            a = a+1;
            build.(fn{f})(:,comp,a) =  sum( diag(loads.(fn{f})(:,comp))*...
                                            X2fr.(fn{f})(:,5:T-8,cond) )';
        end
    end
end

function C = stackStruct(X)
% Will concatenante same fields in diffenret cells
% Cells in X must contain the same fields and dimensions

fStackList = [1:7] ; %Fields number to stack 

ncells = length(X);

if ~isstruct(X)
    C = X{1};
else
    C = X;
end
fn = fieldnames(C); 

for c = 2:ncells
    for f = fStackList
        C.(fn{f}) = vertcat(C.(fn{f}),X{c}.(fn{f}));
    end
end


function X2 = symCells(X)
% Will symmetrises cell by shifting cells to replicated opposite stimuli.

pref  = X.pref;
npref = X.npref;
nbNets = X.nbNets;
 
rmfields = {'Idx','Com','nbTokens','pref','npref','nbNets'};
X        = rmfield(X,rmfields);
 
fn       = fieldnames(X); 
nfields  = length(fn);

N  = size(X.(fn{1}),1);
% T  = size(X.(fn{1}),2);
C  = size(X.(fn{1}),3);

K = N/nbNets; %Nb of neurons per network

% Cloning mirror cells (by shifting cells in each population)
for f = 1:nfields 
    ncirc = length(X.(fn{f}))/K; % number of circshift for this field
    
    for c = 1:C*2
 	% If c is odd, store original data   
        if mod(c,2) 
           X2.(fn{f})(:,:,c) = X.(fn{f})(:,:,(c-1)/2+1);
        
 	% If c is even, shift original data by ~180 degrees orientation
	else
	    % Shifting within each network
            for net = 1:ncirc
                %range of neurons involved 
		range = (net-1)*K+1:K*net; 
		
 		% circshifting network by 180 degres
		X2.(fn{f})(range,:,c) = circshift(X.(fn{f})(range,:,(c-2)/2+1),K/2);
	        net;
 	    end
 
        end
    end
end

