function FR = pcaFormat(nbins,start,stop,Dat,S,Stim,idxStim)
% Function that will put the data into fomart compabtible 
% with PCA algorithm. 
%
%   This include :
%		+ Alignment
%		+ Binning
%		+ concatenation
%		+ Averaging over trials

logi_minCom = abs(Dat.commit)>=start;	% Excludes commitment earlier than 'start' or no commit
logi_right  = Dat.commit>0;		        % Excludes trials where answer is left 

logi = logi_minCom & logi_right; 

dimF = size(Dat.PMD);     % Size of each field
dimC = size(Dat.PMD{1});  % Size of cells in each field

%% Data transformation

% Bin parameters
fn = fields(Dat);		% Field names ; last fields should be fields that don't
                        %               containt neural data

                     
for f = 1:length(fn)-1
    
    goodTr = Dat.(fn{f})(logi);     % Good trials
    
    Dattmp.(fn{f}) = cell(S.nbEx,dimF(2));
    ntT = 1;    %Total number of trial incrementation
    
    FR.(fn{f}) = zeros(dimC(1)*dimF(2),nbins,length(Stim.logIdx));

%% Bin and Alignment 
    for net = 1:dimF(2)
        
        nt = sum(logi(:,net)); %nb of correct trial
        % Removing incorrect trials 
		Dattmp.(fn{f})(1:nt,net) = goodTr(ntT:ntT+nt-1); 
        
        % Only taking relevant sections
        Dattmp.(fn{f})(1:nt,net) = cellfun( @(X) cutTrial(X,start,stop),Dattmp.(fn{f})(1:nt,net), ...
                                    'UniformOutput',false);
        
        % Binning cells
        Dattmp.(fn{f})(1:nt,net) = cellfun( @(X) binmean(X,nbins),Dattmp.(fn{f})(1:nt,net), ...
                                    'UniformOutput',false);        
        ntT = ntT+nt;
    end
    
    
%% Trial classification
    
    if ~S.c | (length(S.c) > 1)
        DatC = classtrial(Dattmp.(fn{f}),Stim,idxStim);
          for c = 1:length(Stim.logIdx)
              maxTrial = max(cellfun('length',DatC(c,:))); % Network with max trial
              DatCtmp  = cell(maxTrial,dimF(2));            
              
              for net = 1:dimF(2)
                  lgt = length(DatC{c,net});
                  DatCtmp(1:lgt,net) = DatC{c,net};
              end
              if ~isempty(DatCtmp)
                 if ~isempty(DatCtmp{1})
                     FR.(fn{f})(:,:,c) = meanTrial(DatCtmp);
                 end
              end
          end
    else 
        FR.(fn{f}) = meanTrial(Dattmp.(fn{f}));
    end
    
end


function D = classtrial(Dat,Stim,idx)
% Will classify trials in respective category if no category specified

dimF = size(Dat);     % Size of each field
dimC = size(Dat{1});  % Size of cells in each field

%Classify stimuli

for c = 1:length(Stim.logIdx)
     trialType{c} = Stim.logIdx{c}(idx);
     for net = 1:dimF(2)
         Dtmp     = Dat(:,net);
         D{c,net} = Dtmp(trialType{c}(:,net)); 
     end
end 
  

function FR = meanTrial(D)

dimF = size(D);     % Size of each field
dimC = size(D{1});  % Size of cells in each field

empt    = cellfun('isempty',D); % idx of empty cells
ntrials = dimF(1)-sum(empt);    % nb of trials per network

A    = cell(1,dimF(2));     % To hold sum of cells
A(:) = {zeros(dimC)};       % Initialize with zeros

D(empt) = {zeros(dimC)};          % Filling empty cells with zeros

for i = 1:dimF(1)
    A = cellfun( @(a,b) a+b, A,D(i,:),...
                 'UniformOutput',false);
end
A  = cellfun(@(a,b) a./b,A,num2cell(ntrials), ...
            'UniformOutput', false); 
FR = cell2mat(A');




		 

