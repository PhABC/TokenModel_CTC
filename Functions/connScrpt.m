function [hnorm,W] =  connScrpt(N,Wsd,Ww,Sunk,spars, wType)
% Creates the tuning curves of population (homogeneous)
r0   = 0;       % Baseline      
rmax = 0.01;    % Peak max
sd   = 30;      % Standart deviation of tuning curve             

%Tuning Curve
hnorm = TuningCurve(r0,rmax,sd,N);

%Connections

% K1 and K2 are internal activity kernel of each region. It is the 
% equivalent of lateral connections within each region. 
W{1,1}  = wMat(Wsd(1,1)*N, Ww(1,1),Sunk(1,1),spars, N, wType);
W{2,2}  = wMat(Wsd(2,2)*N, Ww(2,2),Sunk(2,2),spars, N, wType);

%Weight matrix between regions
W{1,2} = wMat(Wsd(1,2)*N, Ww(1,2), Sunk(1,2),spars, N, wType);
W{2,1} = wMat(Wsd(2,1)*N, Ww(2,1), Sunk(2,1),spars, N, wType);
