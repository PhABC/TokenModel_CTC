function [M1,M2] = CTCsim(S)
%% Simulation

%Unpacking fields
N = S.N; T = S.T; dt = S.dt; tau = S.tau; alpha = S.alpha; beta = S.beta;
gamma = S.gamma; eta = S.eta; Tau = S.Tau; 




    X1 = zeros(N,1);
    X2 = zeros(N,1);
    M1 = zeros(N,T);
    M2 = zeros(N,T);

    Y1 = zeros(N,1);
    Y2 = zeros(N,1);


    for s=1:T/dt
        t = (s-1)*dt;        % Current time in ms

        %activation from PMd1
        s_wY1 = zeros(N,1);
        for j = 1:N
            s_wY1 = s_wY1 + S.w1(:,j).*Y1;
        end

        %activation from PMd2
        s_wY2 = zeros(N,1);
        for j = 1:N
            s_wY2 = s_wY2 + S.w2(:,j).*Y2;
        end

        %within-layer activation 1
        s_KE1 = zeros(N,1);
        for j = 1:N
            s_KE1 = s_KE1 + S.KE(:,j) .* fct(Y1(j));  %nonlinear
    %           s_KE1 = s_KE1 + S.KE(:,j) .* Y1(j);  
        end

        %within-layer activation 2
        s_KE2 = zeros(N,1);
        for j = 1:N
    %         s_KE2 = s_KE2 + S.KE(:,j) .* f(Y2(j)); %nonlinear
              s_KE2 = s_KE2 + S.KE(:,j) .* Y2(j);
        end

        s_I1 = zeros(N,1);
        for j = 1:N
            s_I1 = s_I1 + S.KI(:,j) * fct(Y1(j));
        end

        s_I2 = zeros(N,1);
        for j = 1:N
            s_I2 = s_I2 + S.KI(:,j) * fct(Y2(j));
        end

        Y1 = max(X1-Tau,0);
        Y2 = max(X2-Tau,0);

        E1 = S.V(:,s) + s_wY2 + s_KE1;   
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

end