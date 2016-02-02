function [M1,M2] = CTCsim(S)
%% Simulation

%Unpacking fields
N = S.N; T = S.T; dt = S.dt; tau = S.tau; alpha = S.alpha; beta = S.beta;
gamma = S.gamma; eta = S.eta; Tau = S.Tau; Utype = S.Utype;


% Initialization
    Y1 = zeros(N,N);
    Y2 = zeros(N,N);
    
    X1 = zeros(N,1);
    X2 = zeros(N,1);
    M1 = zeros(N,T);
    M2 = zeros(N,T);


    for s=1:T/dt
        t = (s-1)*dt;        % Current time in ms

        %activation from PMd1
        s_wY1 = sum(S.w1.*Y1,2);

        %activation from PMd2
        s_wY2 = sum(S.w2.*Y2,2);

        %within-layer activation 1    
        s_KE1 = sum(S.KE.*fct(Y1),2); % Nonlinear fct
%       s_KE1 = sum(S.KE.*Y1,2);      % Linear

        %within-layer activation 2
        s_KE2 = sum(S.KE.*fct(Y2),2); % Nonlinear fct
%       s_KE2 = sum(S.KE.*Y2,2);      % Linear
        
        
        s_I1 = sum(S.KI .* fct(Y1),2);
        s_I2 = sum(S.KI .* fct(Y2),2);
        
        %% ???
        
        Y1 = repmat(max(X1-Tau,0),1,N); %% ???
        Y2 = repmat(max(X2-Tau,0),1,N); %% ???

        E1 = s_wY2 + s_KE1 + S.V(:,s);   
        E2 = s_wY1 + s_KE2;

        dX1 = -(alpha.*X1) + (beta - X1).*gamma.*E1 - X1.*s_I1;
        dX2 = -(alpha.*X2) + (beta - X2).*gamma.*E2 - X2.*s_I2;
        
        if Utype == 1 %additive urgency 
            
            dX1 = ( dX1 + S.U(:,s) ) .* tau;
            dX2 = ( dX2 + S.U(:,s) ) .* tau;
        
        elseif Utype == 2 %Multiplicative urgency
        
            dX1 = ( dX1 .* S.U(:,s) ) .* tau;
            dX2 = ( dX2 .* S.U(:,s) ) .* tau;    
        
        end
        
        X1 = X1 + dX1.*dt;
        X2 = X2 + dX2.*dt;
        
        M1(:,s) = X1;
        M2(:,s) = X2;

    end



