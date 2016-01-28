%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD

clear all;
close all;
warning off all;

%% Simulation parameters
S.N     = 50;    % Nb of neurons
S.T     = 1000;  % Simulation time in ms
S.start = 100;   % Start of trial in ms
S.dt    = 1;     % Time step in ms
S.tau   = 0.005; % Time constant

% Stimuli parameters
S.c     = 1;     % type : 1 = easy ~ 2 = misleading ~ 3 = ambiguous
S.nbEx  = 1000;  % Number of stimuli examples to present
S.jumpT = 50;    % interval between each jumps in ms (verify if work with T)
S.stimW = 0.05;  % Amplitude of stimuli

% Input parameters
S.bias = 0.5;   % Additive bias
S.G    = 0.2;   % Noise amplitude

%Expanding time with dt
S.T     = int16(S.T/S.dt);      
S.start = int16(S.start/S.dt);
S.jumpT = int16(S.jumpT/S.dt);

%% Model parameters
S.alpha = 3;     %  
S.beta  = 2;     %
S.gamma = 6;     %
S.eta   = 0.1;   %
S.Tau   = 0.1;   %


%% Tuning curves (homogeneous)
r0   = 0;       %       
rmax = 0.01;    %
sd   = 30;      %                      

S.hnorm = TuningCurve(r0,rmax,sd,S.N);

%% Connections

%Connections matrix
kau = 1.75;
rho = 0.25;
sig = 0.1;

[S.KE, S.KI] = ConMatrix(kau,rho,sig,S.N);

%Weight matrix
[S.w1,S.w2] = WeightMatrix(S.N);


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
S.V = ExtInputs(S,Stim);


%% SIMULATION
[M1,M2] = CTCsim(S);

%find pref units
pref = find(S.V(:,end)>0.5);
npref = setdiff(1:S.N,pref);

%% Figures
% FigureCTC 

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








