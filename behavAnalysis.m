%% Behavioral Analysis 

% Data organisation
Dat.commit    = reshape(commit,S.nbEx*S.nbNet,1);
Dat.comAbs    = abs(Dat.commit);			 %Absolute value of commitment
Dat.token_com = floor( (Dat.comAbs)/S.jumpT+1); %Number of token at commit
idxStim_ok    = FR{1}.Idx;		 %Correct trial idx only

%Correction for nb of token
Dat.token_com(Dat.token_com > 15) = 15;          

for i = 1: S.nbEx
	Dat.stimRawR_com(i,1) = Stim.stimRawR(idxStim(i),Dat.token_com(i)); 
	Dat.probR_com(i,1)    = Stim.probR(idxStim(i),Dat.token_com(i));
end

fn = fields(Dat);

%% Classifying 
for f = 1:length(fn)
    % Correct Trial Classi
    Dat_ok.(fn{f})  = Dat.(fn{f})(Dat.commit>0);	%Correct trial only
    
    % Trials classification
    DatC.(fn{f})    = classTrialStart(Dat.(fn{f})   ,Stim,idxStim);   %All trials
    DatC_ok.(fn{f}) = classTrialStart(Dat_ok.(fn{f}),Stim,idxStim_ok);
end

%% Cumul of success probability
% Success probability of correct trials
cumulProp = cell(1,length(DatC.commit));

for c = 1:length(DatC_ok.commit)
	uni_prop{c} = unique(DatC_ok.probR_com{c});
	nUni = length(uni_prop{c});

	for u = 1:nUni
	    cumulProp{c}(u) = length(DatC_ok.probR_com{c}(DatC_ok.probR_com{c} ... 
                       <= uni_prop{c}(u)))/length(DatC_ok.probR_com{c});		
	end
end
					 

%% Misleading

% Correst misleading vs incorrect misleading 
com_misR =  DatC.commit{2}(DatC.commit{2}>0);	   % Correct answer misleading
com_misL =  abs(DatC.commit{2}(DatC.commit{2}<0)); % Incorrect misleading


%% Plotting 
figure
subplot(2,2,1); hold on
title('Distribution of misleading commits')
hist(abs(DatC.commit{2}),100)
xlabel('Time')

% Reaction time distribution for different class
subplot(2,2,2); hold on
title('Easy vs Misleading vs Ambiguous')
      %Easy
[yeasyH, xeasyH] = hist(DatC.comAbs{1},10);
plot(xeasyH,yeasyH/numel(DatC.commit{1}),'b'); hold on;
      %Misleading
[ymisH,xmisH ]   = hist(abs(DatC.commit{2}),10);
plot(xmisH,ymisH/numel(DatC.comAbs{2}),'r'); hold on;
      %Ambiguous
[yambiH,xambiH] = hist(DatC.comAbs{3},10);
plot(xambiH,yambiH/numel(DatC.comAbs{3}),'g'); hold on;
legend('Easy','Misleading','Ambiguous')
xlabel('Time')

% Misleading correct vs wrong answer
subplot(2,2,3); hold on
title('Correct Misleading vs Incorrect Misleading')
[yMR_H,xMR_H] = hist(com_misR,10); 		
[yML_H,xML_H] = hist(com_misL,10);

plot( xMR_H, yMR_H/numel(yMR_H) ,'g'); hold on; 	   
plot( xML_H, yML_H/numel(yML_H) ,'r'); hold on;	    
legend('Correct','Incorrect')
xlabel('Time')



%Success probability
subplot(2,2,4); hold on;
title('Success prob at decision')
%Easy
plot(uni_prop{1},cumulProp{1}*100,'b'); hold on;
%Misleading
plot(uni_prop{2},cumulProp{2}*100,'r'); hold on;
%Ambiguous
plot(uni_prop{3},cumulProp{3}*100,'g'); hold on;
xlabel('Success probability at decision')
ylabel('Cumulative % of trials')

meanE = round(mean(DatC_ok.probR_com{1})*100);
meanM = round(mean(DatC_ok.probR_com{2})*100);
meanA = round(mean(DatC_ok.probR_com{3})*100);

legend(['Easy       ~ Mean = ', num2str(meanE)],...
       ['Misleading ~ Mean = ', num2str(meanM)], ...
       ['Ambiguous  ~ Mean = ', num2str(meanA)])


%Easy

%[yeasyT,xeasyT] = hist(abs(DatC.stimRawR_com{1}),10);
%plot(xeasyT,yeasyT/numel(DatC.stimRawR_com{1}),'b'); hold on;
%Misleading
%[ymislT,xmislT] = hist(abs(DatC.stimRawR_com{2}),10);
%plot(xmislT,ymislT/numel(DatC.stimRawR_com{2}),'r'); hold on;
%Ambiguous
%[yambiT,xambiT] = hist(abs(DatC.stimRawR_com{3}),10);
%plot(xambiT,yambiT/numel(DatC.stimRawR_com{3}),'g'); hold on;


