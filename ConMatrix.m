function [KE, KI] = ConMatrix(kau,rho,sig,N)
    %% Connections matrix
    
    K = zeros(N);
    
    for i = 1:N
        for j = 1:N

            d = abs(i-j);
            K(i,j) = kau * (exp(-(d^2)/2)/sqrt(2*pi)) - ...
                     0.4 * (exp(-(d^2)/8)/sqrt(2*pi)) - rho;

        end
    end

    KE = max( K,0)+randn(N).*max(max(K)).*0.2;
    KI = max(-K,0)+randn(N).*max(max(K)).*0.2;
    
end

