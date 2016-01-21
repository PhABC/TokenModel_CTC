function [stim,stimInfo] = stimCreation()
%% StimCreation description
%   This function will generate 2^15 unique stimuli based on the token
%   task (see Paul Cisek work). Each stimuli will be composed of 15 tokens
%   which will either to the right or left with an even probability.

%   Stimuli will then be separated in 4 different categories ;

%   Easy:       2T P>60% ~ 5T P>75% ~ 8T P>75% ~ 10T P>80% ~ 13T P>0.9
%   Misleading: 3T P<40% 
%   Ambiguous:  2T P=50% ~ 3T 38%<P<65% ~ 5T 55%<P<65% ~ 7T 55%<P<60%
%   Others:     --------

%   stim     = Matrix containing stimuli (nb Tokens right - nb Tokens left)
%   stimInfo = Structure containing logical indexes and category names 

%% Stimuli creation 
binStim = dec2bin(0:2^15-1);  %All combinations containing 15 T (except 100%)
binStim = binStim - 48;       %Split numbers into different column

count   = triu(ones(15,15));  %to count number of right token at each t
nbRight = binStim*count;      %Nb of tokens to the right at each t
nbLeft  = ~binStim*count;     %Nb of tokens to the left at each t
stim    = nbRight - nbLeft;   %nb Tok right - left at each t


%% Calculating the probability of success
probR          = probRight(nbRight); %Probability if answer right (>0)
stimInfo.probR = probR; 

%% Stimulis classification
%Defining the rules for each category
    %First  row is 'bigger than' conditions
    %Second row is 'equal to' conditions
    %Third  row is 'smaller  than' conditions
 
stimInfo.defName = { 'Easy' 'Misleading' 'Ambiguous' 'Others' } ; 
stimInfo.nbStim   = 2^15;
    
%Easy trials
stimInfo.def{1}  = [ 0, 0.6, 0,0, 0.75, 0,0, 0.75, 0, 0.8, 0,0, 0.9, 0,0  ; ...
                     0,  0,  0,0,   0,  0,0,   0,  0,  0,  0,0,  0,  0,0  ; ...
                     0,  0,  0,0,   0,  0,0,   0,  0,  0,  0,0,  0,  0,0 ];

%Misleading trials                      
stimInfo.def{2}  = [ 0,0,  0,  0,0,0,0,0,0,0,0,0,0,0,0  ; ...
                     0,0,  0,  0,0,0,0,0,0,0,0,0,0,0,0  ; ...
                     0,0, 0.4, 0,0,0,0,0,0,0,0,0,0,0,0 ];

%Ambiguous trials
stimInfo.def{3}  = [ 0,  0, 0.38 ,0, 0.55, 0, 0.55, 0,0,0,0,0,0,0,0  ; ...
                     0, 0.5,  0,  0,   0,  0,   0,  0,0,0,0,0,0,0,0  ; ...
                     0,  0, 0.65, 0, 0.65, 0, 0.60, 0,0,0,0,0,0,0,0 ];

%Classifying                
nonclass = ones(stimInfo.nbStim ,1);

for c = 1:length(stimInfo.def)       %Skipping category 'Others'
    %Conditions for each categories
    B = repmat(stimInfo.def{c}(1,:),stimInfo.nbStim,1); %Bigger than matrix
    B(:,~any(B)) = probR(:,~any(B));  %Take real value when no conditions
    
    E = repmat(stimInfo.def{c}(2,:),stimInfo.nbStim,1); %Equal to matrix
    E(:,~any(E)) = probR(:,~any(E));  %Take real value when no conditions
    
    S = repmat(stimInfo.def{c}(3,:),stimInfo.nbStim,1); %Smaller than matrix
    S(:,~any(S)) = probR(:,~any(S));  %Take real value when no conditions
    
    logMat             = probR >= B & probR == E & probR <= S; %Apply conditions
    stimInfo.logIdx{c} = all(logMat,2);                        %idx with logic array
    stimInfo.idx{c}    = find(stimInfo.logIdx{c} == 1);        %idx with trial nb     
    
    nonclass = nonclass - stimInfo.logIdx{c}; %Every trial that is not classified
end

stimInfo.logIdx{end} = nonclass;              %Saving in the 'Others' class
stimInfo.idx{end}    = find(stimInfo.logIdx{end} == 1);

end

function probR = probRight(Nr)
% Calculates the success probability at each momoment in time
% of selecting the target to the right (assuming right is the correct answer)
%
% Based on following equation :  P(R|NR,NL,NC) = NC! / 2^NC * 
%                                sum_k=0:min(NC,7-NL) ( 1 / (k!(NC-k)!));

Nd     = repmat(1:15,length(Nr),1);        %Nb of token distributed
Nc     = flip(Nd,2)-1;                     %Nb of token in center
Nl     = Nd-Nr;                            %Nb of tokens left

k      = min(Nc,7-Nl);                     %For factorial, nb of sample
k(k<0) = 0;                                %Min boundary,
kSum   = zeros(size(k));                   %Hold the sum of binomial coef


%Calculate the sum of the almost binomial coefficients 
parpool
parfor i = 1:15
       
       subKsum = zeros(length(k),1);      
       for j = 1:length(k) 
           subKsum(j) = sum(1./ (factorial( 0 : k(j,i) ) .* ...
                                 factorial( (Nc(j,i) - (0 : k(j,i) )) )) ); 
       end
       kSum(:,i) = subKsum;
       
end
delete(gcp)

%Calculates the success probabilities
probR   = factorial(Nc) ./ (2.^Nc) .* kSum ;  
      
end

