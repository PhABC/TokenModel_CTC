function [M1,M2,E] = CTCsim(S)
%% Simulation

%Unpacking fields
N = S.N; T = S.T; dt = S.dt; tau = S.tau; alpha = S.alpha; beta = S.beta;
gamma = S.gamma; Tau = S.Tau; Utype = S.Utype; steep = S.steep;


% Initialization
Y1 = zeros(N,1);
Y2 = zeros(N,1);

X1 = zeros(N,1);
X2 = zeros(N,1);
M1 = zeros(N,T);
M2 = zeros(N,T);


for s=1:T/dt
        t = (s-1)*dt;        % Current time in ms

    %activation from PMd1
    s_wY1 = S.W12*fct(Y1,steep_1);

    %activation from PMd2
    s_wY2 = S.W21*fct(Y2,steep_2);

    %Unpacking W matrix into excitation and inhibition
    KE1   =  S.K1;   KE1(KE1<0) = 0;
    KI1   = -S.K1;   KI1(KI1<0) = 0;

    KE2   =  S.K2;   KE2(KE2<0) = 0;
    KI2   = -S.K2;   KI2(KI2<0) = 0;

    %within-layer activation 1    
    s_KE1 = KE1*fct(Y1,steep_1);
    s_KI1 = KI1*fct(Y1,steep_1);

    %within-layer activation 2
    s_KE2 = KE2*fct(Y2,steep_2);
    s_KI2 = KI2*fct(Y2,steep_2);


    %% Calculating activity output and derivatives       
    Y1 = max(X1-Tau,0); 
    Y2 = max(X2-Tau,0); 

    if Utype == 1 %additive urgency 
        
        E1 = s_wY2 + s_KE1 + S.V(:,s) + S.S(:,s) + S.U(:,s);   
        E2 = s_wY1 + s_KE2 + S.V(:,s) + S.U(:,s);

    elseif Utype == 2 %Multiplicative urgency

        E1 = (s_wY2 + s_KE1 + S.V(:,s) + S.S(:,s)).* S.U(:,s);   
        E2 = (s_wY1 + s_KE2 + S.V(:,s)).*S.U(:,s); 

    end
    
    E(:,s) = E1;

    dX1 = -(alpha.*X1) + (beta - X1).*gamma.*E1 - X1.*s_KI1;
    dX2 = -(alpha.*X2) + (beta - X2).*gamma.*E2 - X2.*s_KI2;
    
if s == 300
    s;      %Debug trigger
end
    
        dX1 = dX1 .* tau;
        dX2 = dX2 .* tau;

    X1 = X1 + dX1.*dt;
    X2 = X2 + dX2.*dt;

    M1(:,s) = X1;
    M2(:,s) = X2;

end



