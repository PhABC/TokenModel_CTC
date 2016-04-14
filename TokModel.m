%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMd
    %Each region is refered as R1 and R2.

%% Note : 
% Parallel computing currently not efficient at all. Would need to reduce
% the overhead importantly for it to be interesting. 


%% To do :
% ++++ Optimal parallel computing or GPU computing

%% Global commands
%close all;
clc
warning off all;
clearvars -except seed
tic

% Seed allow you to replay the same trial. Comment out 'seed' if you want to
% reuse the previous trial. 
seed = rng;     %Saving seed (Comment out to reuse previous seed)
rng(seed)         %Loading seed

saveData = 1;    % Will save FR in current dir if == 1

%% Simulation parameters
S.N      = 300;  % Nb of neurons per population
S.T      = 1200; % Simulation time in ms
S.tresh  = 35;	 % Difference between the 2 populations for commitment time

S.nbEx   = 500;   % Number of stimuli examples to present
S.nbNet  = 1;	 % Number of different networks (neurons with different parameters
S.record = 2;	 % if == 1, only involved neurons will be recorded | if == 2, all neurons will be recorded

S.onset  = 120; % onset of trial in ms S.aftcmt = 100; % Stop simulation after X ms S.dt = 1; % Time step in ms S.tau = 0.005; % Time constant
S.aftcmt = 50;  % Time to stop the trial after commitment
S.stRec  = 150; % Time to start recording at. Skipping the first ms of instability in the simulation
S.dt     = 1;    
S.tau    = 0.005; 

S.paralComp = 0;    % run code in parallel if == 1
S.nbWorkers = 2;    % Number of workers for parallel computing
S.pcaForm   = 1;    % Transform data into PCA compabtible format
S.plotting  = 1;    % 1 = plotting ~ 0 = no plots ~ (always plot 1st trial)
S.printDec  = 1;    % Print decision of network
    
S.pcaPlotStSpc = 0;   % Will plot neural state space if True
S.pcaPlotComps = 1:6; % Components to plot in 3d state space
S.pcaPlotConds = 1:6; % Conditions to plot in 3d state space

%% Activation function parameters
%   The steepness of the sigmoid is determined by the parameter 'steep'.
%   Interesting range from 0.04 to 0.64, where 0.4 is near linear regime
%   and 0.64 is a very steep slope. Having a too low value squeez the top
%   and bottom values drasticaly. Could be replaced by pure linear if 
%   real linear regime matters. 

S.RFsteep   = [1,.1]; % 1st dim = if random or not ~ 2nd = standart deviation
S.Fsteep(1) = 0.2 ;   % Steepness of activation function for population 1
S.Fsteep(2) = 0.08  ;  % Steepness of activation function for population 2

S.RFshift   = [1,.01]; % 1st dim = if random or not ~ 2nd = standart deviation
S.Fshift(1) = 60;      % Shifting of the activation function of PMd
S.Fshift(2) = 60;      % Shiftinf of the activation function of M1

%% Connections parameters
%   To have a ff network, make weight from R2 to R1 to 0. .

% Weight between regions

    %Connections R1 to R2 
    S.Ww(1,2)   = .3;   %.3 Amplitude of weight R1 -> R2
    S.Sunk(1,2) = .2;   % Proportion of sunken gaussian. 1 = all inhibitory.   
    S.Wsd(1,2)  = .05;  % 0 < sd < 1 ~ Standart deviation
    
    %Connections R2 to R1
    S.Ww(2,1)   = .15;  %.1
    S.Sunk(2,1) = .2;   %.2
    S.Wsd(2,1)  = .05;  %.05
    
% Weight within regions

    %R1 kernel
    S.Ww(1,1)   = .6;   % Amplitude of weight R1 -> R1 
    S.Sunk(1,1) = .55;  % Proportion of sunken gaussian. 1 = all inhibitory.
    S.Wsd(1,1)  = .1;   % 0 < sd < 1 ~ Standart deviation

    %R2 kernel
    S.Ww(2,2)   = .2;   %.2
    S.Sunk(2,2) = .55;  %.6
    S.Wsd(2,2)  = .1;   %.1
    
%% Input parameters

% Stimuli parameters (defined apriori, WRT stimuli onset, not commit)
S.c     = [1,2,3]; % Classes : 0 = all       | 1 = easy  | 2 = misleading 
                   %           3 = ambiguous | 4 = others   
                   % S.c also accepts vector ( e.g [1.2] )
                    
S.exRan = 1;      % 0 = specified trial ~ 1 = random trial ;
S.exNum = 7528;   % Specifying trial number if non random

S.jumpT = 50;     % interval between each jumps in ms (verify if work with T)
S.stimW = 4.2;     % Amplitude of stimuli ( 0< flip stimuli )

% Bias parameters
S.bias = 10;     % Additive bias strength

% Noise parameters
S.fG   = 10;    % Fast gaussian noise strength (iid)
S.tauFG = 1;    % time constant of fast noise in ms

S.sG   = 0.1;   % Slow gaussian noise strength (shared noise)
S.tauSG = 200;    % time constant of slow noise in ms

% Linear urgency parameters
	%To note, origin and slope will be gaussian distributed for different trials
S.URtrial = [1,.2];
S.URneur  = [1,.2];	% 1st dim = if random or not ~ 2nd = standart deviation

S.Utype  = 1;       % 1 = additive urgency signal ~ 2 = multiplicative urgency signal

S.Uori   = 10;     % origin point for the linear function ~ put 
S.Uslop  = 2.0;   % Slope of the linear urgency function 
S.Uw     = 0.025; % Amplitude of urgency signal [ consider Utype for this value ] 
S.urgmax = 50;

%% Model parameters
S.alpha = 15;      %  Decay factor 
S.beta  = 100;    %  Maximum activity value
S.gamma = 1;      %  Excitation ratio
S.Tau   = 0;      %   

%% Data format parameters for PCA
S.startBin = 400; % Time before commitment 
S.nbins    = 50;  % Number of bins for pca format                       
                  % Time after  commitment is the end of trial

%% Initialization
% Unwrapping certain parameters
N = S.N;  Wsd  = S.Wsd;
Ww= S.Ww; Sunk = S.Sunk;

%Expanding time with dt
S.T      = floor(S.T/S.dt);      
S.onset  = floor(S.onset/S.dt);
S.jumpT  = floor(S.jumpT/S.dt);
S.aftcmt = floor(S.aftcmt/S.dt);
S.stRec  = floor(S.stRec/S.dt);

%Recording and behavioral data 
FR     =  cell(S.nbNet,1);
commit = zeros(S.nbEx,S.nbNet);

%% Stimuli creation
%Creating and saving or loading raw stimuli
if ~exist('Stim','var')
    if exist('StimRaw.mat','file') == 2

        load('StimRaw')

    elseif ~exist('StimRaw.mat','file')

        fprintf('\nCreating and saving stimuli in cd ...\n\n')
        Stim = stimCreation();            % Creating raw stim
        save([pwd '/Data/StimRaw.mat'],'Stim') % Saving variables

    end
end

%% SIMULATION

if S.paralComp
%% Parallel computing
     fprintf('Simulation in parallel ...\n\n')
     
     % Preventing outputs
     S.plotting = 0;   
     S.printDec = 0;   

     S.nbWorkers = min(nbWorkers,S.nbNet);

     if isempty(gcp('nocreate'))
        parpool(nbWorkers);  %open pools
     end

     Gtime = tic;
     parfor net = 1:S.nbNet

         %Creating networks
         %       Creating new networks at every loop with basic
         %      parameters changing (W,tuning, etc.)
         [hnorm,W] =  connScrpt(N,Wsd,Ww,Sunk)
         
         %Stimuli and Bias
         [Uall,stimTrial,idxStim,S] = ExtInputs(S,Stim);

         %Simulation
         fprintf('Simulating %d trials for Network %d ...\n',S.nbEx,net)  
         [FR{net},commit(:,net)] = trialLoop(S,net,hnorm,W,Uall,stimTrial,Stim,idxStim);

     end
     delete(gcp) %close pools
 
else
%% Non-parallel computing    
     fprintf('Simulation ...\n\n')
     Gtime = tic;

    for net = 1:S.nbNet
         if S.nbEx >= 100; S.plotting = 0;  end
         
         %Creating networks
         %      Creating new networks at every loop with basic
         %      parameters changing (W,tuning, etc.)
         [hnorm,W] =  connScrpt(N,Wsd,Ww,Sunk);
         
         %Stimuli and Bias
         [stimTrial,idxStim(:,net),S] = ExtInputs(S,Stim);
         
         
         %Simulation
        fprintf('Simulating %d trials for Network %d of %d ...\n',S.nbEx,net,S.nbNet)  	
        [FR{net},commit(:,net)] = trialLoop(S,net,hnorm,W,stimTrial,Stim,idxStim(:,net),Gtime);

    end
end

commit = (abs(commit)-S.onset).*sign(commit);   %Commit after stim onset 
 
%% PCA
% Plotting neural state space and principal components
if S.pcaPlotStSpc; stateSpace(FR,S.pcaPlotComps,S.pcaPlotConds); end


%% Saving information
if saveData
	Info = S;
	Info.idxStim = idxStim;

	save([pwd '/Data/FR.mat']    ,'FR')
	save([pwd,'/Data/commit.mat'],'commit')
	save([pwd,'/Data/Info.mat']  ,'Info')
end

%% Summary information
fprintf('\nTotal time taken : %s \n',sec2hms(toc(Gtime)))


