function FR = firstJumpFormat(nbins,Dat,Stim,idxStim)

logi_minCom = abs(Dat.commit)>=1;	% Excludes commitment earlier than 'start' or no commit
logi_right  = Dat.commit>0;		% Excludes trials where answer is left 

logi = logi_minCom & logi_right; 

dimC = size(Dat.PMd{1});  % Size of cells in each field

%% Data transformation

% Bin parameters
fn = fields(Dat);		% Field names ; last fields should be fields that don't
                        %               containt neural data

for f = 1:length(fn)-2
    
    goodTr  = Dat.(fn{f})(logi);     % Good trials

    FR.(fn{f}) = zeros(dimC(1),nbins,length(Stim.logIdx));

    %% Bin and Alignment 

    % adding NaN at end of trials
    goodTr = nanfill(goodTr);
     
    % Binning cells
    goodTr  = cellfun( @(X) binmean(X,nbins),goodTr, ...
                       'UniformOutput',false);
                   
    goodIdx  =  idxStim(logi);
    goodCom  =  Dat.commit(logi);
    goodNTok =  Dat.nbTokens(logi);
  
%% Trial classification
      
   %DatC = classTrialCommit(goodTr, Stim, goodIdx,goodNTok);  %wrt to commit
    DatC =  classTrialStart(goodTr, Stim, goodIdx);           %wrt onset

for c = 1:length(Stim.logIdx)              
       if ~isempty(DatC{c})
          if ~isempty(DatC{1})
              FR.(fn{f})(:,:,c) = meanTrial(DatC{c});                 
          end
       end
    end   
end



function FR = meanTrial(D)

% Mean of all trials in cell
stackD  = cat(3,D{:});
FR      = nanmean(stackD,3);



function FR = nanfill(D)

maxLen = max(cell2mat(cellfun(@length, D, 'uni',false)));

for c = 1:length(D)
    D{c}(:,length(D{c}):maxLen) = NaN; 
end
FR = D; 
