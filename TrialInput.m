function [ stim,V,U,pref,npref ] = TrialInput(S,trial)
%% TrialInput
% Will put all inputs in right format for each trial in the simulation   

%Unpacking fields
  T=S.T; N = S.N; nbEx=S.nbEx; onset=S.onset; stimW=S.stimW; Utype=S.Utype;

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
  hstim  = S.hnorm .* repmat(tokG,1,N);             
  pref   = sum(hstim)>0;
  npref  = sum(hstim)<0;
  
%Cummulative input for each neurons based on stim and pref 
  stim   = hstim'   * repmat(S.stim(trial,:),360,1);

%Negative values indicate their non prefered direction
  stim(stim<0) = 0;

  
%% Fast noise
  FG = randn(N,T);

%% Bias, fast noise and slow noise
  V  = S.bias + repmat(S.SG(trial,:),N,1)*S.sG + FG*S.fG;
  U  = repmat(S.urg(trial,:),N,1);
  
  
  
 
  
  
  
  
  
  
  
  
  
  
   
   
   
   
   