%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD

% clear all;
close all;
warning off all;

 seed = rng;       %Saving seed
%rng(seed)           %Loading seed

%% Simulation parameters
S.N      = 50;    % Nb of neurons
S.T      = 1000;  % Simulation time in ms
S.onset  = 100;   % onset of trial in ms
S.dt     = 1;     % Time step in ms
S.tau    = 0.005; % Time constant

%% Input parameters

% Stimuli parameters
S.c      = 1;     % type : 1 = easy ~ 2 = misleading ~ 3 = ambiguous
S.nbEx   = 1;  % Number of stimuli examples to present
S.jumpT  = 50;    % interval between each jumps in ms (verify if work with T)
S.stimW  = 0.01;  % Amplitude of stimuli

% Bias parameters
S.bias   = 1;   % Additive bias strength

% Noise parameters
S.fG     = 0.01;  % Fast gaussian noise strength (iid)
S.sG     = 0.2 ;  % Slow gaussian noise strength (shared noise)

% Linear urgency parameters
	%To note, origin and slope will be gaussian distributed for different trials
S.Utype = 1;	% 1 = additive urgency signal ~ 2 = multiplicative urgency signal
S.Uori  = 1;    % origin point for the linear function ~ put 
S.Uslop = 1;    % Slope of the linear urgency function 
S.Uw    = 1;	% Amplitude of urgency signal [ consider Utype for this value ] 


%% Model parameters
S.alpha = 3;     %  
S.beta  = 2;     %
S.gamma = 6;     %
S.eta   = 0.1;   %
S.Tau   = 0.1;   %

%Expanding time with dt
S.T     = floor(S.T/S.dt);      
S.onset = floor(S.onset/S.dt);
S.jumpT = floor(S.jumpT/S.dt);

%% Tuning curves (homogeneous)
r0   = 0;       % Baseline      
rmax = 0.01;    % Peak max
sd   = 30;      % Standart deviation of tuning curve             

S.hnorm = TuningCurve(r0,rmax,sd,S.N);



%% Connections

%Connections matrix
kau = 1.75;
rho = 0.25;
sig = 0.1;

[S.KE, S.KI] = ConMatrix(kau,rho,sig,S.N);

%Weight matrix
[S.w1,S.w2]  = WeightMatrix(S.N);


%% Stimuli and Bias

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

PD = 270;
OD = 90;
sd = 30;

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








