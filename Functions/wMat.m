function W = wMat(sd,Ww,sunkR,N)
    
    x = 1:N*2;
    all_shifts = linspace(1,N*2,N*2);    
    all_g = zeros(N*2,N);
    k = 1;

    for i = 1:N
        s = all_shifts(k);
        all_g(:,k) = exp((-(x-s).^2)./(2*sd.^2));
        k = k+1;
    end
    
    % Excitation connections
    E = all_g(1:N,:)+all_g(N+1:end,:)+... 
           flip(flip(all_g(N+1:end,:),2));
    
    % Cummulative weight matrix
    W = E*Ww;
    
    sunk = max(max(W))*sunkR;
       
    W = W-sunk + randn(size(W))*Ww*0.2;    
    W = W.*100/N;			%Normalize wrt 100 as parameters were found with 100 N
end
