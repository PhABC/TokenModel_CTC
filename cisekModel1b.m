%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD
    %Each region is refered as R1 and R2.

% ++++ SD for inhibition not implemented yet ?
% ++++ proper output format for PCA decision manifold 
% ++++ Clip urgency so it doesn't go crazy at the end

% close all;
warning off all;
clearvars -except seed

% Seed allow you to replay the same trial. Comment out 'seed' if you want to
% reuse the previous trial. 
  seed = rng;       %Saving seed (Comment out to reuse previous seed)
  rng(seed)         %Loading seed

%% Simulation parameters
S.N      = 200;  % Nb of neurons per population
S.T      = 1000;  % Simulation time in ms
S.onset  = 120;   % onset of trial in ms
S.dt     = 1;     % Time step in ms
S.tau    = 0.005; % Time constant
plotting = 1;     % 1 = plotting ~ 0 = no plots

%% Activation function parameters
%   The steepness of the sigmoid is determined by the parameter 'steep'.
%   Interesting range from 0.04 to 0.64, where 0.4 is near linear regime
%   and 0.64 is a very steep slope. Having a too low value squeez the top
%   and bottom values drasticaly. Could be replaced by pure linear if 
%   real linear regime matters. 
 
S.steep_1 = 0.16; %Steepness of activation function for population 1
S.steep_2 = 0.04; %Steepness of activation function for population 2

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
    S.Ww_12   = .01;   % Amplitude of weight R1 -> R2
    S.Sunk_12 =  0.2;  % Proportion of sunken gaussian. 1 = all inhibitory.   
    S.Wsd_12  = .05;    % 0 < sd < 1 ~ Standart deviation
    
    %Connections R2 to R1
    S.Ww_21   =  .0;   
    S.Sunk_21 =  .2;   
    S.Wsd_21  =  .05;   
    
% Weight within regions

    %R1 kernel
    S.Ww_11   =  .5;   % Amplitude of weight R1 -> R1 
    S.Sunk_11 =  .6;   % Proportion of sunken gaussian. 1 = all inhibitory.
    S.Wsd_11  =  .1;   % 0 < sd < 1 ~ Standart deviation

    %R2 kernel
    S.Ww_22   = .15;   
    S.Sunk_22 = .6; 
    S.Wsd_22  = .1;   
    

%% Input parameters

% Stimuli parameters
S.c      = 2;    % type : 1 = easy ~ 2 = misleading ~ 3 = ambiguous
S.nbEx   = 10;    % Number of stimuli examples to present
S.jumpT  = 50;   % interval between each jumps in ms (verify if work with T)
S.stimW  = 3;    % Amplitude of stimuli ( 0< flip stimuli )

% Bias parameters
S.bias   = 5;   % Additive bias strength

% Noise parameters
S.fG     = 10;   % Fast gaussian noise strength (iid)
S.sG     = 0.1;  % Slow gaussian noise strength (shared noise)

% Linear urgency parameters
	%To note, origin and slope will be gaussian distributed for different trials
S.Urand  = 0;	 % 1 = random slope every trial ~ 0 = same slope every trial
S.Utype  = 1;	 % 1 = additive urgency signal ~ 2 = multiplicative urgency signal
S.Uori   = 100;  % origin point for the linear function ~ put 
S.Uslop  = 3;     % Slope of the linear urgency function 
S.Uw     = 0.015;    % Amplitude of urgency signal [ consider Utype for this value ] 


%% Model parameters
S.alpha = 10;    %  Decay factor 
S.beta  = 100;   %  Maximum activity value
S.gamma = 1;     %  Excitation ratio
S.Tau   = 0;     %  

% Unwrapping certain parameters
N = S.N;


%% Tuning curves (homogeneous)
r0   = 0;       % Baseline      
rmax = 0.01;    % Peak max
sd   = 30;      % Standart deviation of tuning curve             

S.hnorm = TuningCurve(r0,rmax,sd,N);


%% Connections
%Connections matrix

% K1 and K2 are internal activity kernel of each region. It is the 
% equivalent of lateral connections within each region. 
S.K1  = wMat(S.Wsd_11*N, S.Ww_11,S.Sunk_11, N);
S.K2  = wMat(S.Wsd_22*N, S.Ww_22,S.Sunk_22, N);

%Weight matrix between regions
S.W12 = wMat(S.Wsd_12*N, S.Ww_12, S.Sunk_12, N);
S.W21 = wMat(S.Wsd_21*N, S.Ww_21, S.Sunk_21, N);

%% Stimuli and Bias

%Expanding time with dt
S.T     = floor(S.T/S.dt);      
S.onset = floor(S.onset/S.dt);
S.jumpT = floor(S.jumpT/S.dt);

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

%Creating inputs
fprintf('\nExternal input creation ... \n');
[ S.urg,S.stim,S.SG ] = ExtInputs(S,Stim);

%% SIMULATION
fprintf('Simulation ...\n\n')

for trial = 1:S.nbEx
    
    fprintf('Trial %d of %d \n',trial,S.nbEx) 
    % stimuli preferences are also assigned in the following fct
    [ S.S, S.V, S.U, S.pref, S.npref ] = TrialInput(S,trial);
    [ M1,M2 ] = CTCsim(S);
    
%   Figures
    if plotting == 1 
	FigureCTC
    end

end


 

return;

%PCA
% [Up Sp Vp] = pca(M1(pref,:));   %Was pca2 before 
% [Un Sn Vn] = pca(M1(npref,:));  %Was pca2 before 
% figure;
% h = plot3(Vp(:,1),Vp(:,2),Vp(:,3),'Color',[0 128 0]./255);
% set(h,'linewidth',8);
% hold on;
% h = plot3(Vn(:,1),Vn(:,2),Vn(:,3),'Color',[128 0 0]./255);
% set(h,'linewidth',8);
% set(gca,'FontSize',24);
% xlabel('PCA1');
% ylabel('PCA2');
% zlabel('PCA3');
% 
% figure;
% [Up Sp Vp] = pca2(M2(pref,:));
% [Un Sn Vn] = pca2(M2(npref,:));
% h = plot3(Vp(:,1),Vp(:,2),Vp(:,3),'Color',[128 255 128]./255);
% set(h,'linewidth',8);
% hold on;
% h = plot3(Vn(:,1),Vn(:,2),Vn(:,3),'Color',[255 128 128]./255);
% set(h,'linewidth',8);
% set(gca,'FontSize',24);








