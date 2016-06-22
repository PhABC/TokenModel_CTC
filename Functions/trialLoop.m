function [FR,allcommit]  = trialLoop(S,net,hnorm,W,stimTrial,Stim,idxStim,Gtime)
% Simulation with all trials

tic
batchSize = 1000; % Number of trials per batch ~ to save RAM

%Will split trials to save RAM
nBatch  = floor(S.nbEx/batchSize);  %Number of batches of size batchSize
resi    = mod(S.nbEx,batchSize);    %Residual trials for last batch

% Defining all batches size
if resi > 0
    batches = horzcat(ones(1,nBatch)*batchSize,resi);
else
    batches = ones(1,nBatch)*batchSize;             
end

nBatch   = length(batches);
totTrial = 1; %Total number of trials

% Will contain all commits and concatenate them
allcommit = [];

NetState{1,1} = zeros(S.N,1); NetState{1,2} = zeros(S.N,1);
NetState{2,1} = zeros(S.N,1); NetState{2,2} = zeros(S.N,1);

for batch = 1:nBatch
    %Number of trial in the batch
    ntrial = batches(batch); 
    
    %Initializing
    Dat.PMd = cell(ntrial,1);
    Dat.M1  = cell(ntrial,1);
    commit   = zeros(ntrial,1);
    nbTokens = zeros(ntrial,1);

    idxStimBatch = idxStim(totTrial:totTrial+ntrial-1);
    
    for trial = 1:ntrial
        totTrial     = batchSize*(batch-1) + trial; 
            
        % stimuli preferences are also assigned in the following fct
        [ St,V,U,pref,npref ] = TrialInput(S,trial,hnorm,stimTrial);
        [ PMd_,M1_,commit(trial),St,NetState] = CTCsim(S,W,St,V,U,pref,npref,NetState);     

        if S.printDec 

            %Number of tokens (15 is the max)
            nbTokens(trial) =  floor((abs(commit(trial)/S.dt)-S.onset)/S.jumpT)+1;
            if nbTokens(trial) > 15; nbTokens(trial) = 15; end
            
            fprintf('      Trial :   %d of %d    ~    Decision :   ',totTrial,S.nbEx)

            if commit(trial) > 0
                fprintf('Correct   (tokens: %d)\n', nbTokens(trial));
            elseif commit(trial) < 0
                fprintf('Incorrect (tokens: %d)\n', nbTokens(trial));
            else 
                fprintf('--- \n')
            end

        end

        %Figures
        if S.plotting || trial == 1 && batch == 1  
           FigureCTC
        end

        if S.nbEx > 100
            if ~mod(totTrial,S.nbEx/10)
                fprintf('\n       -------------------------------------------------------\n')
                fprintf('             Net #%d                Time : %s\n', net, sec2hms(toc(Gtime)));
                fprintf('       -------------------------------------------------------\n\n')
             end
        end

        %Window of time to record
        timeCom = floor(S.stRec:min(abs(commit(trial)/S.dt)+S.aftcmt,S.T));

        %Sample recording
        if S.record == 1
        %Only outputting values of prefered neurons
            Dat.PMd{trial} = PMd_(logical(pref+npref),timeCom);
            Dat.M1{trial}  =  M1_(logical(pref+npref),timeCom); 
        elseif S.record == 2
        %Outputting values of all neurons
            Dat.PMd{trial} = PMd_(:,timeCom);
            Dat.M1{trial}  =  M1_(:,timeCom);
        end

%         if S.pcaForm == 1
%               FR = pcaFormat(S.nbins,S.startBin,S.aftcmt,Dat,Stim,idxStim);
%         else
%               FR = Dat;
%         end  
    end

%     Dat.commit   = vertcat(Dat.commit,commit);
%     Dat.nbTokens = vertcat(Dat.nbTokens,nbTokens);

    Dat.commit   = commit;
    Dat.nbTokens = nbTokens;

    if S.pcaForm == 1
        FR{batch} = pcaFormat(S.nbins,S.startBin,S.aftcmt,Dat,Stim,idxStimBatch);
    else
        FR{batch} = Dat;
    end
    allcommit = vertcat(allcommit,commit);
    
    FR{batch}.pref  = pref';
    FR{batch}.npref = npref';
end

% Batch averaging weighted by number of trials per batch
if S.pcaForm && nBatch > 1
	
	ntList = zeros(nBatch,1); %nb of trials for each batch

	for b = 1:nBatch
		ntList(b)=length(FR{b}.Idx);
	end

	ntTot  = sum(ntList); %Total number of good trials	
	
	%Weighting average and stacking
	FRtmp = FR{1};
	FRtmp.PMd = FR{1}.PMd*(ntList(1)/ntTot);
	FRtmp.M1  = FR{1}.M1 *(ntList(1)/ntTot);
	
	for b = 2:nBatch
	    FRtmp.PMd(:,:,:,b) = FR{b}.PMd*(ntList(b)/ntTot);
	    FRtmp.M1(:,:,:,b)  = FR{b}.M1*(ntList(b)/ntTot);
	    FRtmp.Idx = vertcat(FRtmp.Idx,FR{b}.Idx);
	    FRtmp.Com = vertcat(FRtmp.Com,FR{b}.Com);
	    FRtmp.nbTokens = vertcat(FRtmp.nbTokens,FR{b}.nbTokens);
	end
	
	FRtmp.PMd = nansum(FRtmp.PMd,4); % Summation part of weighted averaging
	FRtmp.M1 = nansum(FRtmp.M1,4);	 % Summation part of weighted averaging

	FR = FRtmp;
end



function [A1,A2,commit,St,NetState] = CTCsim(S,W,St,V,U,pref,npref,NetState) %% Simulation
%% Running the simulation

%Unpacking fields 
N = S.N; T = S.T; dt = S.dt; tau = S.tau; alpha = S.alpha; beta = S.beta;
gamma = S.gamma; Tau = S.Tau; Utype = S.Utype; tresh=S.tresh;

% Initialization 

% If initializing wiht previous trial network state
if S.contiTr 
    X1 = NetState{1,1}; X2 = NetState{1,2};
    Y1 = NetState{2,1}; Y2 = NetState{2,2};
else
    X1 = zeros(N,1); X2 = zeros(N,1);
    Y1 = zeros(N,1); Y2 = zeros(N,1);
end

A1     = zeros(N,T); A2 = zeros(N,T);
commit = 0; %If commit stays 0, it means the network didn't cross tresh

for s=1:T
    t = (s-1)*dt; % Current time in ms

    PMd_A  = fct(Y1,S.STEEP(:,1),S.SHIFT(:,1));  % PMd activity after transfer function
    M1_A   = fct(Y2,S.STEEP(:,2),S.SHIFT(:,2));  % M1  activity after transfer function
    
    %activation from PMd1 to M1
    s_wY1 = W{1,2}*PMd_A;

    %activation from M1 to PMd
    s_wY2 = W{2,1}*M1_A ;

    %Unpacking W matrix into excitation and inhibition
    %   Could be out of the loop with diff name
    KE1 =  W{1,1}; KE1(KE1<0) = 0;
    KI1 = -W{1,1}; KI1(KI1<0) = 0;

    KE2 =  W{2,2}; KE2(KE2<0) = 0;
    KI2 = -W{2,2}; KI2(KI2<0) = 0;

    %within-layer activation 1
    s_KE1 = KE1*PMd_A;
    s_KI1 = KI1*PMd_A;

    %within-layer activation 2
    s_KE2 = KE2*M1_A;
    s_KI2 = KI2*M1_A ;

    %% Calculating activity output and derivatives
    Y1 = max(X1-Tau,0);
    Y2 = max(X2-Tau,0);

    if Utype == 1 %additive urgency
        
        E1 = s_wY2 + s_KE1 + V(:,s) + St(:,s) + U.PMd(:,s);
        E2 = s_wY1 + s_KE2 + V(:,s) + U.M1(:,s)  + St(:,s);

    elseif Utype == 2 %Multiplicative urgency

        E1 = ( s_wY2 + s_KE1 + St(:,s) ).*U.PMd(:,s)+ V(:,s) ;
        E2 = ( s_wY1 + s_KE2 ).*U.M1(:,s) + V(:,s) ;

    end
    
    dX1 = -(alpha.*X1) + (beta - X1).*gamma.*E1 - X1.*s_KI1;
    dX2 = -(alpha.*X2) + (beta - X2).*gamma.*E2 - X2.*s_KI2;

if t == 500;
1;
end %Debug trigger
    
    dX1 = dX1 .* tau;
    dX2 = dX2 .* tau;

    X1 = X1 + dX1*dt;
    X2 = X2 + dX2*dt;

 % Min value of 0
    X1(X1<0) = 0;
    X2(X2<0) = 0;

    A1(:,s) = X1;
    A2(:,s) = X2;

%% Treshold
    % Used to use mean, but mean is about 8x slower than max and
    % qualitatively results are identical.
    diffPop = mean(X1(pref)) - mean(X1(npref));

    if abs(diffPop) >= tresh && ~commit
        %Decaying inputs to network 
       if S.IDecay; [U,St] = decayInput(S,U,St,s); end
       commit = sign(sign(S.stimW)*diffPop)*(t); %Sign indicate direction
    elseif s == T-S.aftcmt 
        %Decaying inputs to network
       if S.IDecay; [U,St] = decayInput(S,U,St,s); end      
    end

    %To let run the simulation a bit longer.
    if s == abs(floor(commit/S.dt))+S.aftcmt && commit ~= 0 || s == T
       %Saving end state of network
       NetState{1,1} = X1; NetState{1,2} = X2;
       NetState{2,1} = Y1; NetState{2,2} = Y2;
	    break
    end

end


function [U,St] = decayInput(S,U,St,s)
% Will decay inputs to birng the network back
% to the initial equilibrium

   decayT = S.aftcmt-50; %Time for inputs to decay

   %Decay of urgency for PMd
   downU_PMd = linspaceMat(U.PMd(:,s),U.PMd(:,1), decayT);
   U.PMd(:,s:s+decayT-1) = downU_PMd; 
   U.PMd(:,s+decayT:end) = U.PMd(1,1);

   %Decay of urgency for PMd   
   downU_M1 = linspaceMat(U.M1(:,s),U.M1(:,1), decayT);
   U.M1(:,s:s+decayT-1) = downU_M1; 
   U.M1(:,s+decayT:end) = U.M1(1,1);
  
   %Decay of stimuli input
   downS = linspaceMat(St(:,s),St(:,1),decayT);
   St(:,s:s+decayT-1) = downS;
   St(:,s+decayT:end) = St(1,1);

