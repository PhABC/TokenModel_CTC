function D = classTrialStart(Dat,Stim,idx)
% Will classify trials in respective category if no category specified

%Classify stimuli
for c = 1:length(Stim.logIdx)
     trialType{c} = Stim.logIdx{c}(idx);
             D{c} = Dat(trialType{c});   
end 