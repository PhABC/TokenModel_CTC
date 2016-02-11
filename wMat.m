function W = wMat(r0,rmax,sd,WEw,WIw,N)
    
    x = 1:N*2;
    all_shifts = linspace(1,N*2,N*2);    
    all_g = zeros(N*2,N);
    k = 1;

    for i = 1:N
        s = all_shifts(k);    
        all_g(:,k) = r0 + (rmax.*exp((-(x-s).^2)./(2*sd.^2)));
        k = k+1;
    end
    
    % Excitation connections
    E = all_g(1:N,:)+all_g(N+1:end,:)+... 
           flip(flip(all_g(N+1:end,:),2));
    
    % Inhibition connections
    I = circshift(E,N/2,2);
    
    % Cummulative weight matrix
    W = E*WEw - I*WIw;
   
    
    W = W + randn(size(W))*WEw*.05;    

end
