function tuningMat = TuningCurve(r0,rmax,sd,N)
%% TUNING CURVE
    % Generate a tuning curve based on gaussian
    
    x = 1:720;
    all_shifts = linspace(1,720,N*2);    
    all_g = zeros(720,N);
    k = 1;

    for i = 1:N
        s = all_shifts(k);    
        all_g(:,k) = r0 + (rmax.*exp((-(x-s).^2)./(2*sd.^2)));
        k = k+1;
    end
    
    tuningMat = all_g(1:360,:)+all_g(361:end,:);
    
end
