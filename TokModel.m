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
% seed = rng;     % Saving seed (Comment out to reuse previous seed)
rng(seed)         % Loading seed

saveData = 1;     % Will save FR in current dir if == 1

%% Simulation parameters
S.N      = 300;  % Nb of neurons per layer
S.T      = 1000; % Simulation time in ms
S.tresh  = 35;	 % Difference between the 2 populations for commitment time

S.nbEx    = 1000;  % Number of stimuli examples to present
S.nbNet   = 1; 	  % Number of different networks (neurons with different parameters
S.record  = 2;	  % if == 1, only core neurons will be recorded | if == 2, all neurons will be recorded
S.pcaForm = 1;    % Transform data into PCA compabtible format

S.onset  = 100; % onset of trial in ms
S.aftcmt = 50;  % Time to stop the trial after commitment
S.stRec  = 100; % Time to start recording at. Skipping the first ms of instability in the simulation
S.dt     = 1;    
S.tau    = 0.005; 

S.plotting  = 1; % 1 = plotting ~ 0 = no plots ~ (always plot 1st trial)
S.printDec  = 1; % Print decision of network
S.paralComp = 0; % UNDER DEVELOPMENT. LEAVE AT 0. ~ Run code in parallel if == 1
S.nbWorkers = 2; % Number of workers for parallel computing
    
S.pcaPlotStSpc = 0;   % Will plot neural state space if True
S.pcaPlotComps = 1:6; % Components to plot in 3d state space
S.pcaPlotConds = 1:6; % Conditions to plot in 3d state space

%% Activation function parameters
%   The steepness of the sigmoid is determined by the parameter 'steep'.
%   Interesting range from 0.04 to 0.64, where 0.4 is near linear regime
%   and 0.64 is a very steep slope. Having a too low value squeez the top
%   and bottom values drasticaly. Could be replaced by pure linear if 
%   real linear regime matters. 

S.RFsteep   = [1,.02]; % 1st dim = if random or not ~ 2nd = standart deviation
S.Fsteep(1) = 0.12 ;   % Steepness of activation function for population 1
S.Fsteep(2) = 0.12 ;   % Steepness of activation function for population 2

S.RFshift   = [1,2];   % 1st dim = if random or not ~ 2nd = standart deviation
S.Fshift(1) = 60;      % Shifting of the activation function of PMd
S.Fshift(2) = 60;      % Shiftinf of the activation function of M1

%% Connections parameters
% Note:
%       + Connectivity is a function of tuning curve similarity with
%              added noise.
%       + Distribution parameters for each type can be found in the
%              script Functions/wMat.m

% Connectivity sparsness proportion
S.spars = 0.7; % ( 1 == "fully connected" )
	       
% Connectivity distribution types
%	0 : Sunken gaussian distribution
% 	1 : Opposite tuning inhibition
%	2 : Mexicain hat distribution (i.e. near inhibition)
% 	3 : Near and opposite tuning inhibition 
S.wType = 0;

% Weight between regions
%        To have a ff network, make weight from R2 to R1 to 0

    %Connections R1 to R2 
    S.Ww(1,2)   = .2;   % Amplitude of weight R1 -> R2
    S.Sunk(1,2) = .2;   % Proportion of sunken gaussian. 1 = all inhibitory.   
    S.Wsd(1,2)  = .1;   % 0 < sd < 1 ~ Standart deviation
    
    %Connections R2 to R1
    S.Ww(2,1)   = .2;
    S.Sunk(2,1) = .2;   
    S.Wsd(2,1)  = .1;  
    
% Weight within regions

    %R1 kernel
    S.Ww(1,1)   = .5;   % Amplitude of weight R1 -> R1 
    S.Sunk(1,1) = .55;  % Proportion of sunken gaussian. 1 = all inhibitory.
    S.Wsd(1,1)  = .1;   % 0 < sd < 1 ~ Standart deviation

    %R2 kernel
    S.Ww(2,2)   = .5;   
    S.Sunk(2,2) = .55;  
    S.Wsd(2,2)  = .1;   
    
%% Input parameters

% Stimuli parameters (defined apriori, WRT stimuli onset, not commit)
S.c     = [1,2,3]; % Classes : 0 = All       | 1 = Easy   | 2 = Misleading 
                   %           3 = Ambiguous | 4 = Others   
                   % S.c also accepts vector ( e.g [1.2] )
                    
S.exRan = 1;       % 0 = specified trial ~ 1 = random trial
S.exNum = 7528;    % Specifying trial number if non random

S.jumpT = 50;      % interval between each jumps in ms (verify if work with T)
S.stimW = 4;       % Amplitude of stimuli ( sign flips stimuli )

% Bias parameters
S.bias  = 2;    % Additive bias amplitude

% Noise parameters
S.fG    = 10;   % Fast gaussian noise strength (iid)
S.tauFG = 2;    % time constant of fast noise in ms

S.sG    = 0.1;  % Slow gaussian noise strength (shared noise)
S.tauSG = 200;  % time constant of slow noise in ms (has to be a multiple of T)

% Linear urgency parameters
S.Utuning = 1;          % 1: urgency as a function of stimuli tuning
                        % To change the coefficient mean or biais, see
                        % S.USdist variable below. More details is given 
                        % as to how to build the urgency distribution.

S.URtrial = [1,.2];     % Random mean urgency between trials
S.URneur  = [1,.4];     % Random urgency signal between neurons
	                    % [X,~] = random or not (1 or 0) ~ [~,X] = absolute standart deviation

S.Utype  = 1;     % 1 = additive urgency signal ~ 2 = multiplicative urgency signal
        		  % One needs to change the parameters value and test for stability 
		          % if switching to multiplicative urgency signal.


S.Uori   = 0;     % origin point for the linear function ~ put 
S.Uslop  = 4.5;   % Slope of the linear urgency function 
S.Uw     = 0.01;  % Amplitude of urgency signal [ consider Utype for this value ] 
S.urgmax = 60;

%% Model parameters
S.alpha = 15;     %  Decay factor 
S.beta  = 100;    %  Maximum activity value
S.gamma = 1;      %  Excitation ratio
S.Tau   = 0;      %   

%% Data format parameters for PCA
S.startBin = 350; % Time before commitment 
S.nbins    = 50;  % Number of bins for pca format                       
                  % Time after  commitment is the end of trial

%% Initialization
% Unwrapping certain parameters
N = S.N;  Wsd  = S.Wsd;  spars = S.spars;
Ww= S.Ww; Sunk = S.Sunk;

%Expanding time with dt
S.T      = floor(S.T/S.dt);      
S.onset  = floor(S.onset/S.dt);
S.jumpT  = floor(S.jumpT/S.dt);
S.aftcmt = floor(S.aftcmt/S.dt);
S.stRec  = floor(S.stRec/S.dt);

%Recording and behavioral data 
FR     = cell(S.nbNet,1);
commit = zeros(S.nbEx,S.nbNet);

%% Stimuli creation
%Creating and saving or loading raw stimuli
if ~exist('Stim','var')
    if exist('StimRaw.mat','file') == 2

        load('StimRaw')

    elseif ~exist('StimRaw.mat','file')

        fprintf('\nCreating and saving stimuli in cd ...\n\n')
        Stim = stimCreation();            % Creating raw stim
        save([pwd '/TokenModel_CTC/Data/StimRaw.mat'],'Stim') % Saving variables

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
         [hnorm,W] =  connScrpt(N,Wsd,Ww,Sunk,spars,S.wType)
         
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
         [hnorm,W] = connScrpt(N,Wsd,Ww,Sunk,spars,S.wType);
     
         % Urgency tuning    
	     if S.URneur(1)  
            % First argument of randParam is the mean of tuning coefficient for each region.
	    	%      0 < USdist implies more positive tuning
	        %      Note: +It seems important to have negatively tuned urgency in order to have a realistic
            %                   distribution of eigenvalues for the components found by PCA.
            %            +If S.Utuning is 1, S.Udist will be added to another matrix of coefficient
            %            +If S.Utuning is 0, recommended value is between 0.2 and 0.5
            %            +If S.Utuning is 1, recommended value is around -0.5, as it is added value
            %                   going up to 2.5.            
            S.USdist  = randParam([-0.5,-0.5],S.URneur(2),N); 	
	 	    
            % Standardizing urgency tuning by mean
	 	    S.USdist  = S.USdist./repmat(max(abs(S.USdist)),N,1); 
	     else
		    S.USdist = repmat([0.2,0.2],N,1);	    
	     end
             
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
if S.pcaPlotStSpc; stateSpace(FR,0,S.pcaPlotConds,S.pcaPlotComps); end


%% Saving information
if saveData
	Info = S;
	Info.idxStim = idxStim;

	save([pwd '/TokenModel_CTC/Data/FR.mat']    ,'FR')
	save([pwd,'/TokenModel_CTC/Data/commit.mat'],'commit')
	save([pwd,'/TokenModel_CTC/Data/Info.mat']  ,'Info')
end

%% Summary information
fprintf('\nTotal time taken : %s \n',sec2hms(toc(Gtime)))


