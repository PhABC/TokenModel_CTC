function [ urg,stim,SG ] = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli
%Unpacking fields
   T=S.T; dt=S.dt; nbEx=S.nbEx; c=S.c;jumpT=S.jumpT; 
   Uslop=S.Uslop; onset=S.onset; Uw=S.Uw; Uori=S.Uori;

     
%Initializing
   nstim   = Stim.nbStim{c};               %Number of stim of chosen c
   stimRaw = Stim.stimRawR(Stim.idx{c},:); %Trials foc condition c
   
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
   
%Choosing examples
   idx     = randi(nstim,[nbEx,1]); %Select nbEx random integer within nstim
   stim    = stimAll(idx,:);
   
%% Urgency   
%Baseline urgency signal
   Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
   urg  = zeros(nbEx,T);  
   urg(:,onset:end) = repmat(linspace(0,Uend,T-onset+1),... 
                             [nbEx,1]);      % Urgency starts at stim onset
%    urg  = urg/max(urg(:,end));                             % Normalized


%Creating gauss distribution   
if S.Urand
   urg  = urg .* abs(repmat(randn(nbEx,1)+Uslop,1,T))...  % Random slopes
               + repmat(0.1*randn(nbEx,1)+Uori,1,T);      % Random origins
end
%    urg  = urg/max(urg(:,end));                            % Normalized
   urg(urg<0) = 0;                                        % Del values < 0

   urg = urg*Uw;

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
