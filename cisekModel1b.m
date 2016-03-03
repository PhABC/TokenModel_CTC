%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD
    %Each region is refered as R1 and R2.

%% Note : 
% Script is currently formated to ease writing, but this increases 
% the stress on the ram for paralell computing.proportionnaly to the nb of
% neurons. This can be change with a little effort.

%% To do :
% ++++ SD for inhibition not implemented yet ?
% ++++ proper output format for PCA decision manifold 

%% Global commands
% close all;
clc
warning off all;
clearvars -except seed
tic

% Seed allow you to replay the same trial. Comment out 'seed' if you want to
% reuse the previous trial. 
  seed = rng;       %Saving seed (Comment out to reuse previous seed)
  rng(seed)         %Loading seed

%% Simulation parameters
S.N      = 200;   % Nb of neurons per population
S.T      = 2000;  % Simulation time in ms
S.tresh  = 35;	  % Difference between the 2 populations for commitment time

S.nbEx   = 100;   % Number of stimuli examples to present
S.nbNet  = 2;	  % Number of different networks (neurons with different parameters)
    
S.onset  = 120;   % onset of trial in ms
S.aftcmt = 100;   % Stop simulation after X ms
S.dt     = 1;     % Time step in ms
S.tau    = 0.005; % Time constant

paralComp  = 0;   % run code in parallel if == 1
nbWorkers  = 4;   % Number of workers for parallel computing
saveFR     = 1;   % Will save FR in current dir if == 1
S.plotting = 1;   % 1 = plotting ~ 0 = no plots
S.printDec = 1;   % Print decision of network

%% Activation function parameters
%   The steepness of the sigmoid is determined by the parameter 'steep'.
%   Interesting range from 0.04 to 0.64, where 0.4 is near linear regime
%   and 0.64 is a very steep slope. Having a too low value squeez the top
%   and bottom values drasticaly. Could be replaced by pure linear if 
%   real linear regime matters. 
 
S.steep(1) = 0.16; %Steepness of activation function for population 1
S.steep(2) = 0.04; %Steepness of activation function for population 2

%% Connections parameters
%   To have a ff network, make weight from R2 to R1 to 0. Can adjust the
%   strength of excitation and inhibition by changing the amplitude ( Ww )
%   and the proportion of the guassian that is negative (sunk).
%
%   The standar devations (SD) have to be between 0 and 1 as they are
%   proportional to the number of neurons. More specifically, this sd value
%   is multiplied by the number of neurons.

% Weight between regions

    %Connections R1 to R2 
    S.Ww(1,2)   = .3;   % Amplitude of weight R1 -> R2
    S.Sunk(1,2) = .2;   % Proportion of sunken gaussian. 1 = all inhibitory.   
    S.Wsd(1,2)  = .05;  % 0 < sd < 1 ~ Standart deviation
    
    %Connections R2 to R1
    S.Ww(2,1)   =  .0;   
    S.Sunk(2,1) =  .2;   
    S.Wsd(2,1)  =  .05;   
    
% Weight within regions

    %R1 kernel
    S.Ww(1,1)   =  .5;   % Amplitude of weight R1 -> R1 
    S.Sunk(1,1) =  .6;   % Proportion of sunken gaussian. 1 = all inhibitory.
    S.Wsd(1,1)  =  .1;   % 0 < sd < 1 ~ Standart deviation

    %R2 kernel
    S.Ww(2,2)   = .15;   
    S.Sunk(2,2) = .6; 
    S.Wsd(2,2)  = .1;   
    
%% Input parameters

% Stimuli parameters
S.c      = 2;    % type : 0 = all | 1 = easy | 2 = misleading | 3 = ambiguous | 4 = others                  
S.exRan  = 1;	 % 0 = specified trial ~ 1 = random trial ;
S.exNum  = 330;  % Specifyin trial number if non random

S.jumpT  = 50;   % interval between each jumps in ms (verify if work with T)
S.stimW  = 4;    % Amplitude of stimuli ( 0< flip stimuli )

% Tuning curves (homogeneous)
r0   = 0;       % Baseline      
rmax = 0.01;    % Peak max
sd   = 30;      % Standart deviation of tuning curve             


% Bias parameters
S.bias   = 5;    % Additive bias strength

% Noise parameters
S.fG     = 10;   % Fast gaussian noise strength (iid)
S.sG     = 0.1;  % Slow gaussian noise strength (shared noise)

% Linear urgency parameters
	%To note, origin and slope will be gaussian distributed for different trials
S.Urand  = 0;	 % 1 = random slope every trial ~ 0 = same slope every trial
S.Utype  = 1;	 % 1 = additive urgency signal ~ 2 = multiplicative urgency signal

S.Uori   = 0;    % origin point for the linear function ~ put 
S.Uslop  = 3;    % Slope of the linear urgency function 
S.Uw     = 0.015;    % Amplitude of urgency signal [ consider Utype for this value ] 
S.urgmax = 30;

%% Model parameters
S.alpha = 15;    %  Decay factor 
S.beta  = 100;   %  Maximum activity value
S.gamma = 1;     %  Excitation ratio
S.Tau   = 0;     %  

%% Data format parameters
startBin = 300;	 % Time before commitment 
nbins    = 80;                        
% Time after  commitment is the end of trial

%% Initialization
% Unwrapping certain parameters
N = S.N; Wsd=S.Wsd;
Ww=S.Ww; Sunk=S.Sunk;

%Expanding time with dt
S.T     = floor(S.T/S.dt);      
S.onset = floor(S.onset/S.dt);
S.jumpT = floor(S.jumpT/S.dt);

%Creating matrices for general results
% Dat.full   =  cell(S.nbNet,S.nbEx);
M1     =  cell(S.nbEx,S.nbNet);
PMD    =  cell(S.nbEx,S.nbNet);
commit = zeros(S.nbEx,S.nbNet);

%% Stimuli creation
%Creating and saving or loading raw stimuli
if ~exist('Stim','var')
    if exist('StimRaw.mat','file') == 2

        load('StimRaw')

    elseif ~exist('StimRaw.mat','file')

        fprintf('\nCreating and saving stimuli in cd ...\n\n')
        Stim = stimCreation();            % Creating raw stim
        save([pwd '/StimRaw.mat'],'Stim') % Saving variables

    end
end

%% Creating networks
% Creating new networks at every loop with basic parameters changing (W,tuning, etc.)

fprintf('Creating %d networks ... \n\n',S.nbNet);
for net = 1:S.nbNet

    %% Tuning Curve
 	S.hnorm{net} = TuningCurve(r0,rmax,sd,N);

	%% Connections
	%Connections matrix

	% K1 and K2 are internal activity kernel of each region. It is the 
	% equivalent of lateral connections within each region. 
 	S.W{net}{1,1}  = wMat(Wsd(1,1)*N, Ww(1,1),Sunk(1,1), N);
 	S.W{net}{2,2}  = wMat(Wsd(2,2)*N, Ww(2,2),Sunk(2,2), N);

	%Weight matrix between regions
 	S.W{net}{1,2} = wMat(Wsd(1,2)*N, Ww(1,2), Sunk(1,2), N);
 	S.W{net}{2,1} = wMat(Wsd(2,1)*N, Ww(2,1), Sunk(2,1), N);
% 
	%% Stimuli and Bias

	%Creating inputs
 	[ S.urg{net},S.stim{net},S.SG{net},S.idxStim(:,net)] = ExtInputs(S,Stim);
end


	%% SIMULATION
	
 
 if paralComp
     fprintf('Simulation in parallel ...\n\n')
     % Preventing outputs
     S.plotting = 0;   
     S.printDec = 0;   
     
     nbWorkers = min(nbWorkers,S.nbNet);
     
     parpool(nbWorkers);  %open pools
     parfor net = 1:S.nbNet

        fprintf('\nSimulating %d trials for Network %d ...\n',S.nbEx,net)  
        [PMD(:,net),M1(:,net),commit(:,net)] = trialLoop(S,net);

     end
     delete(gcp) %close pools
    
 else
     fprintf('Simulation ...\n\n')
     for net = 1:S.nbNet
        
      fprintf('\nSimulating %d trials for Network %d ...\n',S.nbEx,net)  
          [PMD(:,net),M1(:,net),commit(:,net)] = trialLoop(S,net);

     end
     
 end
 
commit     = (abs(commit)-S.onset).*sign(commit);   %Commit after stim onset 
 
Dat.PMD    = PMD;
Dat.M1     =  M1;
Dat.commit = commit;

FR = pcaFormat(nbins,startBin,S.aftcmt,Dat,S,Stim);

%% Aligning trials with respect to commitment
%FR = pcaFormat(nbins,startBin,S.aftcmt,Dat,S) 	    
	    
 
fprintf('\nTotal time taken in seconds : %f \n',toc)


if saveFR
save([pwd '/FR.mat'],'FR')
end




