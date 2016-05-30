function [ stim,V,U,pref,npref ] = TrialInput(S,trial,hnorm,stimAll)
%% TrialInput
% Will put all inputs in right format for each trial in the simulation   

%Unpacking fields
T=S.T; N = S.N; stimW=S.stimW; dt = S.dt; Uw=S.Uw; Uori = S.Uori;
USdist = S.USdist; urgmax=S.urgmax; Uslop=S.Uslop;onset = S.onset;

%% Stimuli and orientation preference
%Orientation for token 
sd     = 30;      % Standart deviation of token orientation
m      = 180;     % Mean for center of direction 

%Because stimuli are opposite directions, the left is represented with
%negative values and the right stim are positive. If values have same
%sign as the actual stimuli, then the stim value will be positive,
%otherwise it will be negative. Removing negative values later.

gaussL = -Gauss2d([360,1],[0,m],sd,1);  % 1D Distribution - Left angle    
gaussR = -circshift(gaussL,-m);         % Translation for right angle

tokG   = gaussL + gaussR;               % Combining left & right  

%Neurons preference with token orientation
htune  = hnorm .* repmat(tokG,1,N);            
pref   = sum(htune)> .4;
npref  = sum(htune)<-.4;

%Cummulative input for each neurons based on stim and pref 
stim   = htune'   * repmat(stimAll(trial,:),360,1) * stimW;

%Negative values indicate their non prefered direction
stim(stim<0) = 0;

%% Fast noise
noiseFG = randn(N,T/(S.tauFG/dt)+1);

stepFG = S.tauFG/dt; 

FG = zeros(N,T);
i = 1; 
%Expanding slow noise
for t = 1:stepFG:T
    FG(:,t:t+stepFG-1) = linspaceMat(noiseFG(:,i),noiseFG(:,i+1),stepFG);
    i=i+1;
end

%% Slow noise

noiseSG = randn(1,T/(S.tauSG/dt)+1);

stepSG = S.tauSG/dt; 

SG = zeros(1,T);
i = 1; 
%Expanding slow noise
for t = 1:stepSG:T
    SG(1,t:t+stepSG-1) = linspaceMat(noiseSG(i),noiseSG(i+1),stepSG);
    i=i+1;
end

%     %Smoothing out noise
%     for ex = 1:nbEx
%         SG(ex,:) = smooth(SG(ex,:),100/dt);
%     end

%% Bias, fast noise and slow noise
V  = S.bias + repmat(SG,N,1)*S.sG + FG*S.fG;
   
%% Urgency   
%Baseline urgency signal
Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
urg  = zeros(1,T);  
urg(1:onset)   = Uori;
urg(onset:end) = repmat(linspace(Uori,Uend,T-onset+1),... 
                         [1,1]);      % Urgency starts at stim onset

%Creating gauss distribution   
if S.URtrial(1)
   urg  = urg  + urg  .* repmat(S.URtrial(2)*randn,1,T) + ...  % Random slopes
          Uori + Uori .* repmat(S.URtrial(2)*randn,1,T);       % Random origins
end
   
Uall.PMd = cell(1,1);
Uall.M1 = cell(1,1);

% Applying urgency affinity  
% The urgency signal tuning can:
% 	+Random 
% 	+Function of stimuli tuning solely
% 	+Function of sitmuli tuning with randomness
% Randonmess is set by S.URneur(1) in the TokModel script.
% Urgency tuned wrt stimuli prefererence is set by 

%       Add abs() if only positive urgency signal


% Urgency tuning as a function of stimuli tuning
if S.Utuning
    %Tuning affinity based on stimuli preferences of neurons
    Utuning = ( sum(abs(htune))/max(sum(abs(htune)))...  %Normalizing as a function of stim preference
              + mean(sum(abs(htune))) )';

    % Random urgency tuning
    Uall.PMd = repmat(urg,N,1) .* repmat(USdist(:,1)+Utuning,1,T);
     Uall.M1 = repmat(urg,N,1) .* repmat(USdist(:,2)+Utuning,1,T); 
else
    % Random urgency tuning
    Uall.PMd = repmat(urg,N,1) .* repmat(USdist(:,1),1,T);
     Uall.M1 = repmat(urg,N,1) .* repmat(USdist(:,2),1,T); 
end


%Urgency weigth
Uall = structfun( @(X) X*Uw, Uall, 'UniformOutput',false);

%Setting boundary 
% Max value
Uall.PMd(Uall.PMd >  urgmax) =  urgmax; %max urg 
 Uall.M1(Uall.M1  >  urgmax) =  urgmax; %max urg 
% Min value
Uall.PMd(Uall.PMd < -urgmax) = -urgmax; %max urg 
 Uall.M1(Uall.M1  < -urgmax) = -urgmax; %max urg 

% To only take positive urgency signals
% Uall = structfun( @(X) abs(X), Uall, 'UniformOutput',false);
 
%Output structure
U.PMd = Uall.PMd;
U.M1  = Uall.M1;

 
  
  
  
  
  
  
  
  
  
  
   
   
   
   
   
