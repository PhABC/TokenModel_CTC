function [ stim,V,U,pref,npref] = TrialInput(S,trial,hnorm,stimAll)
%% TrialInput
% Will put all inputs in right format for each trial in the simulation   

%Unpacking fields
T=S.T; N = S.N; stimW=S.stimW; dt = S.dt; Uw=S.Uw; Uori = S.Uori;
urgmax=S.urgmax; Uslop=S.Uslop;onset = S.onset;

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
hstim  = hnorm .* repmat(tokG,1,N);             
pref   = sum(hstim)> .4;
npref  = sum(hstim)<-.4;

%Cummulative input for each neurons based on stim and pref 
stim   = hstim'   * repmat(stimAll(trial,:),360,1) * stimW;

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
          Uori + Uori .* repmat(S.URtrial(2)*randn,1,T);    	             % Random origins
end

urg(urg<0) = 0;                                        % Del values < 0

urg = urg*Uw;
urg(urg>urgmax) = urgmax; %max urg 

   
Uall.PMd = cell(1,1);
Uall.M1 = cell(1,1);

if S.URneur(1) == 1
    
   % PMd & M1 slopes factors
   sl = randParam([1,1],S.URneur(2),N); 
    
    % PMD & M1 slopes
    Uall.PMd = repmat(urg,N,1) .* repmat(sl(:,1),1,T);
     Uall.M1 = repmat(urg,N,1) .* repmat(sl(:,2),1,T); 
    
end 

U.PMd = Uall.PMd;
U.M1  = Uall.M1;


  
 
  
  
  
  
  
  
  
  
  
  
   
   
   
   
   
