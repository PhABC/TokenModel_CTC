function [ urg,stim,SG,idx ] = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli
%Unpacking fields
   T=S.T; dt=S.dt; nbEx=S.nbEx; c=S.c;jumpT=S.jumpT; 
   Uslop=S.Uslop; onset=S.onset; Uw=S.Uw;
   urgmax=S.urgmax; exRan=S.exRan; exNum=S.exNum;
     
%Initializing
if ~c
   nstim   = Stim.nbStimR;               %Number of stim of chosen c
   stimRaw = Stim.stimRawR; %Trials for condition c
elseif length(c) ==1   
   nstim   = Stim.nbStim{c};               %Number of stim of chosen c
   stimRaw = Stim.stimRawR(Stim.idx{c},:); %Trials for condition c
elseif length(c) > 1
    nstim   = Stim.nbStim{c(1)};
    stimRaw = Stim.stimRawR(Stim.idx{c(1)},:); 
    for cc = c(2:end)
        nstim   = nstim+Stim.nbStim{cc};
        stimRaw = vertcat(stimRaw,Stim.stimRawR(Stim.idx{cc},:));
    end
end
   
%Logic matrix to expand stimuli wrt T
   timeLogiMat = zeros(15,T);
   i = 1;
   
   for bin = onset:jumpT:T
       if i < 15
         timeLogiMat(i,bin:bin+jumpT-1) = 1;
       else
         timeLogiMat(i,bin:end) = 1;
         break;
       end
       i = i+1;
   end
   
   stimAll = stimRaw*timeLogiMat;   %Stimuli from c expanded wrt T
   
   
%% IDX IS NOT THE IDX WITH RESPECT TO ALL STIM EXCEPT IF C = 0 (i.e. using all stim)  
%Choosing examples
   if ~exRan
        idx = exNum;
   else
        idx = randi(nstim,[nbEx,1]); %Select nbEx random integer within nstim
   end
   stim = stimAll(idx,:);
   
%% Urgency   
%Baseline urgency signal
   Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
   urg  = zeros(nbEx,T);  
   urg(:,onset:end) = repmat(linspace(0,Uend,T-onset+1),... 
                             [nbEx,1]);      % Urgency starts at stim onset
%    urg  = urg/max(urg(:,end));                             % Normalized


%Creating gauss distribution   
if S.Urand
   urg  = urg  + urg .* abs(repmat(randn(nbEx,1)*.5,1,T)) + ...  % Random slopes
          uori + repmat(0.1*randn(nbEx,1),1,T);      % Random origins
end
%    urg  = urg/max(urg(:,end));                            % Normalized
   urg(urg<0) = 0;                                        % Del values < 0

   urg = urg*Uw;
   urg(urg>urgmax) = urgmax; %max urg 

%% Slow noise
    tauSG = 500;                          % time constant of slow noise in ms
    noise = randn(nbEx,T/(tauSG/dt)+1);
    
    step = tauSG/dt; 
    
    SG = zeros(nbEx,T);
    i = 1; 
    
    %Expanding slow noise
    for t = 1:step:T
        SG(:,t:t+step-1) = linspaceMat(noise(:,i),noise(:,i+1),step);
        i=i+1;
    end
    
    %Smoothing out noise
    for ex = 1:nbEx
        SG(ex,:) = smooth(SG(ex,:),100/dt);
    end
    
end
