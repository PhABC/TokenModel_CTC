%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD
    %Each region is refered as R1 and R2.

% ++++ SD for inhibition not implemented yet. 
% ++++ change transpose operation in CTCsim (slow down simulation) 

    
% clear all;
% close all;
warning off all;
clearvars -except seed
 
  seed = rng;       %Saving seed
 rng(seed)           %Loading seed

%% Simulation parameters
S.N      = 100;    % Nb of neurons
S.T      = 1000;  % Simulation time in ms
S.onset  = 120;   % onset of trial in ms
S.dt     = 1;     % Time step in ms
S.tau    = 0.005; % Time constant

%% Connections parameters
%   To have a ff network, make weight from R2 to R1 to 0. Can adjust the
%   strength of excitation and inhibition with WEw and WIw respectively.
%
%   The standar devations (SD) have to be between 0 and 1 as they are
%   proportional to the number of neurons. More specifically, this sd value
%   is multiplied by the number of neurons. **NOT IMPLEMENTED FOR
%   INHIBITION**

% Weight between regions

    %Connections R1 to R2 
    S.WEw_12  = .02;   % Amplitude of excitory   weight R1 -> R2
    S.WIw_12  = .000;  % Amplitude of inhibitory weight R1 -> R2
    
    S.WEsd_12 = .1;   % 0 < sd < 1 ~ Standart deviation of E
%   S.WIsd_12 = .1;   % 0 < sd < 1 ~ Standart deviation of I

    %Connections R2 to R1
    S.WEw_21  = 0.0;   % Amplitude of excitory   weight R2 -> R1
    S.WIw_21  = 0.0;   % Amplitude of inhibitory weight R2 -> R1
    
    S.WEsd_21 = .1;   % 0 < sd < 1 ~ Standart deviation of E
%   S.WIsd_21 = .1;   % 0 < sd < 1 ~ Standart deviation of I

% Weight within regions

    %R1 kernel
    S.WEw_1   = .01;   % Amplitude of excitory   weight R1 -> R1
    S.WIw_1   = .1;   % Amplitude of inhibitory weight R1 -> R1
    
    S.WEsd_1  = .1;   % 0 < sd < 1 ~ Standart deviation of E
%   S.WIsd_1  = .1;   % 0 < sd < 1 ~ Standart deviation of I

    %R2 kernel
    S.WEw_2   = .001;   % Amplitude of excitory   weight R2 -> R2
    S.WIw_2   = .04;   % Amplitude of inhibitory weight R2 -> R2
    
    S.WEsd_2  = .1;   % 0 < sd < 1 ~ Standart deviation of E
%   S.WIsd_2  = .1;   % 0 < sd < 1 ~ Standart deviation of I

%% Input parameters

% Stimuli parameters
S.c      = 2;    % type : 1 = easy ~ 2 = misleading ~ 3 = ambiguous
S.nbEx   = 1;    % Number of stimuli examples to present
S.jumpT  = 50;   % interval between each jumps in ms (verify if work with T)
S.stimW  = 5;  % Amplitude of stimuli ( 0< flip stimuli )

% Bias parameters
S.bias   = 10;    % Additive bias strength

% Noise parameters
S.fG     = 20;   % Fast gaussian noise strength (iid)
S.sG     = 0.1;  % Slow gaussian noise strength (shared noise)

% Linear urgency parameters
	%To note, origin and slope will be gaussian distributed for different trials
S.Utype  = 1;	 % 1 = additive urgency signal ~ 2 = multiplicative urgency signal
S.Uori   = 1;    % origin point for the linear function ~ put 
S.Uslop  = 1;    % Slope of the linear urgency function 
S.Uw     = 20;  % Amplitude of urgency signal [ consider Utype for this value ] 


%% Model parameters
S.alpha = 5;     %  
S.beta  = 100;     %  %Maximum activity value
S.gamma = 1;     %
S.eta   = 0;   %
S.Tau   = 0;   %

% Unwrapping certain parameters
N = S.N;


%% Tuning curves (homogeneous)
r0   = 0;       % Baseline      
rmax = 0.01;    % Peak max
sd   = 30;      % Standart deviation of tuning curve             

S.hnorm = TuningCurve(r0,rmax,sd,N);


%% Connections
%Connections matrix
% kau = 1.75;
% rho = 0.25;
% sig = 0.1;
% [S.KE, S.KI]    = ConMatrix(kau,rho,sig,N);

% K1 and K2 are internal activity kernel of each region. It is the 
% equivalent of lateral connections within each region. 
S.K1  = wMat(0,1, S.WEsd_1*N, S.WEw_1, S.WIw_1, N);
S.K2  = wMat(0,1, S.WEsd_2*N, S.WEw_2, S.WIw_2, N);

%Weight matrix between regions
S.W12 = wMat(0,1, S.WEsd_12*N, S.WEw_12, S.WIw_12, N);
S.W21 = wMat(0,1, S.WEsd_21*N, S.WEw_21, S.WIw_21, N);



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
    
end

%% Figures
FigureCTC 

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








