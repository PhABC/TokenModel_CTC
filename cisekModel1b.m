%% SIMPLIFIED VERSION WITH ONLY a 2-layer PMD

clear all;
close all;
warning off all;

%% Simulation parameters
N = 50;         % Nb of neurons
T = 500;        % Simulation time
dt = 1;         % Time step in milliseconds
tau = 0.005;    % 


%% Model parameters
alpha = 3;      %  
beta = 2;       %
gamma = 6;      %
eta = 0.1;      %
Tau = 0.1;      %


%% Tuning curves (homogeneous)
r0 = 0;
rmax = 0.01;
sd = 30;                             %was 10
x = 1:720;
all_shifts = linspace(1,720,N*2);    
all_g = zeros(720,N);
k = 1;

for i = 1:N
    s = all_shifts(k);    
    all_g(:,k) = r0 + (rmax.*exp((-(x-s).^2)./(2*sd.^2)));
    k = k+1;
end
hnorm = all_g(1:360,:)+all_g(361:end,:);

%% Connections generation

%Weight matrix
w = zeros(N);
for i = 1:N
    w(i,max(i-3,1):i) = linspace(0,0.4,length(max(i-3,1):i));
    w(i,i:min(i+3,N)) = linspace(0.4,0,length(i:min(i+3,N)));
end
w1 = w+randn(N).*max(max(w)).*0.01;
% w2 = w+randn(N).*max(max(w)).*0.01; %recurrent connections
w2 = zeros(N); %FF connections
 
%Connections matrix
kau = 1.75;
rho = 0.25;
sig = 0.1;
K = zeros(N);
for i = 1:N
    for j = 1:N
        d = abs(i-j);
        K(i,j) = kau * (exp(-(d^2)/2)/sqrt(2*pi)) - ...
            0.4 * (exp(-(d^2)/8)/sqrt(2*pi)) - rho;
    end
end

KE = max(K,0)+randn(N).*max(max(K)).*0.2;
KI = max(-K,0)+randn(N).*max(max(K)).*0.2;

%HERE: ADD SOME INPUT (V)
V = zeros(N,T);

%% Stimuli and Bias
%easy, non-misleading stimulus
PD = 270;
OD = 90;
stim = zeros(360,T);
stim(PD,1:T) = 1:T;
x = 1:360;
for t = 1:T
    stim(:,t) =  (t/T).*exp((-(x-PD).^2)./(2*sd.^2)) - ...
        (t/T).*exp((-(x-OD).^2)./(2*sd.^2));
    k = k+1;
end
V = hnorm' * stim;
V = V+repmat([1:T]./1000,N,1);  %urgency
% imagesc(V)
% return;

%misleading stimulus
% PD = 270;
% OD = 90;
% stim = zeros(360,T);
% stim(PD,1:T) = 1:T;
% x = 1:360;
% for t = 1:T-100
%     stim(:,t) =  (t/T).*exp((-(x-OD).^2)./(2*sd.^2)) - ...
%         (t/T).*exp((-(x-PD).^2)./(2*sd.^2));
%     k = k+1;
% end
% 
% for t = T-100:T
%     stim(:,t) =  (t/T).*exp((-(x-PD).^2)./(2*sd.^2)) - ...
%         (t/T).*exp((-(x-OD).^2)./(2*sd.^2));
%     k = k+1;
% end
% V = hnorm' *stim;
% imagesc(V)
% return;

X1 = zeros(N,1);
X2 = zeros(N,1);
M1 = zeros(N,T);
M2 = zeros(N,T);

Y1 = zeros(N,1);
Y2 = zeros(N,1);

%MAIN FOR LOOP
for s=1:T/dt
    t = (s-1)*dt;        % Current time in ms
    
    %activation from PMd1
    s_wY1 = zeros(N,1);
    for j = 1:N
        s_wY1 = s_wY1 + w1(:,j).*Y1;
    end
    
    %activation from PMd2
    s_wY2 = zeros(N,1);
    for j = 1:N
        s_wY2 = s_wY2 + w2(:,j).*Y2;
    end
    
    %within-layer activation 1
    s_KE1 = zeros(N,1);
    for j = 1:N
        s_KE1 = s_KE1 + KE(:,j) .* fct(Y1(j));  %nonlinear
%           s_KE1 = s_KE1 + KE(:,j) .* Y1(j);  
    end
    
    %within-layer activation 2
    s_KE2 = zeros(N,1);
    for j = 1:N
%         s_KE2 = s_KE2 + KE(:,j) .* f(Y2(j)); %nonlinear
          s_KE2 = s_KE2 + KE(:,j) .* Y2(j);
    end
    
    s_I1 = zeros(N,1);
    for j = 1:N
        s_I1 = s_I1 + KI(:,j) * fct(Y1(j));
    end
    
    s_I2 = zeros(N,1);
    for j = 1:N
        s_I2 = s_I2 + KI(:,j) * fct(Y2(j));
    end
    
    Y1 = max(X1-Tau,0);
    Y2 = max(X2-Tau,0);
    
    E1 = V(:,s) + s_wY2 + s_KE1;   
    E2 = s_wY1 + s_KE2;
    
    dX1 = -(alpha.*X1) + (beta - X1).*gamma.*E1 - X1.*s_I1+(eta.*randn(N,1));
    dX2 = -(alpha.*X2) + (beta - X2).*gamma.*E2 - X2.*s_I2+(eta.*randn(N,1));
    dX1 = dX1.*tau;
    dX2 = dX2.*tau;
    X1 = X1 + dX1.*dt;
    X2 = X2 + dX2.*dt;
    M1(:,s) = X1;
    M2(:,s) = X2;

end

%find pref units
pref = find(V(:,end)>0.5);
npref = setdiff(1:N,pref);

figure;
hold on;
plot(mean(M1(npref,:)),'r','linewidth',4);
plot(mean(M1(pref,:)),'g','linewidth',4);
set(gca,'Linewidth',3);
set(gca,'FontSize',40);
xlabel('time (ms)');
ylabel('activation');
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








