function [A1,A2,commit] = CTCsim(S,net,St,V,U,pref,npref) %% Simulation

%Unpacking fields 
N = S.N; T = S.T; dt = S.dt; tau = S.tau; alpha = S.alpha; beta = S.beta;
gamma = S.gamma; Tau = S.Tau; Utype = S.Utype; steep = S.steep; tresh=S.tresh;

% Initialization 
Y1 = zeros(N,1); Y2 = zeros(N,1);
X1 = zeros(N,1); X2 = zeros(N,1);
A1 = zeros(N,T); A2 = zeros(N,T);

commit = 0; %If commit stays 0, it means the network didn't cross tresh

for s=1:T
        t = (s-1)*dt; % Current time in ms

    PMD_A  = fct(Y1,steep(1));  % PMD activity after transfer function
    M1_A   = fct(Y2,steep(2));  % M1  activity after transfer function
    
    %activation from PMd1 to M1
    s_wY1 = S.W{net}{1,2}*PMD_A;

    %activation from M1 to PMd
    s_wY2 = S.W{net}{2,1}*M1_A ;

    %Unpacking W matrix into excitation and inhibition
    KE1 =  S.W{net}{1,1}; KE1(KE1<0) = 0;
    KI1 = -S.W{net}{1,1}; KI1(KI1<0) = 0;

    KE2 =  S.W{net}{2,2}; KE2(KE2<0) = 0;
    KI2 = -S.W{net}{2,2}; KI2(KI2<0) = 0;

    %within-layer activation 1
    s_KE1 = KE1*PMD_A;
    s_KI1 = KI1*PMD_A;

    %within-layer activation 2
    s_KE2 = KE2*M1_A;
    s_KI2 = KI2*M1_A ;

    %% Calculating activity output and derivatives
    Y1 = max(X1-Tau,0);
    Y2 = max(X2-Tau,0);

    if Utype == 1 %additive urgency
        
        E1 = s_wY2 + s_KE1 + V(:,s) + St(:,s) + U(:,s);
        E2 = s_wY1 + s_KE2 + V(:,s) + U(:,s);

    elseif Utype == 2 %Multiplicative urgency

        E1 = ( s_wY2 + s_KE1 + V(:,s) + St(:,s) ).*U(:,s);
        E2 = ( s_wY1 + s_KE2 + V(:,s) ).*U(:,s);

    end
    
    dX1 = -(alpha.*X1) + (beta - X1).*gamma.*E1 - X1.*s_KI1;
    dX2 = -(alpha.*X2) + (beta - X2).*gamma.*E2 - X2.*s_KI2;

if t == 500;
    s; 
end %Debug trigger
    
    dX1 = dX1 .* tau;
    dX2 = dX2 .* tau;

    X1 = X1 + dX1*dt;
    X2 = X2 + dX2*dt;

    A1(:,s) = X1;
    A2(:,s) = X2;

%% Treshold
    % Used to use mean, but mean is about 8x slower than max and
    % qualitatively results are identical.
    diffPop = max(X1(pref)) - max(X2(npref));

    if abs(diffPop) >= tresh && ~commit
	commit = sign(diffPop)*(t); %Sign indicate direction	
    end

    %To let run the simulation a bit longer.
    if t == abs(commit)+S.aftcmt && commit ~= 0
	    break
    end

end



