%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD

clear all;
close all;
warning off all;

%% Simulation parameters

S.N   = 50;     % Nb of neurons
S.T   = 500;    % Simulation time
S.dt  = 1;      % Time step in milliseconds
S.tau = 0.005;  % 


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
PD = 270;
OD = 90;
sd = 30;

S.V = ExtInputs(PD,OD,sd,S.hnorm,S.T,S.N);


%% SIMULATION
[M1,M2] = CTCsim(S);

%find pref units
pref = find(S.V(:,end)>0.5);
npref = setdiff(1:S.N,pref);

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








