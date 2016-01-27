function V = ExtInputs(PD,OD,sd,hnorm,T,N)
%% External inputs
    %All the inputs not coming from the neural activity
    
%easy, non-misleading stimulus
stim = zeros(360,T);
stim(PD,1:T) = 1:T;
x = 1:360;
k=1;
for t = 1:T
    stim(:,t) =  (t/T).*exp((-(x-PD).^2)./(2*sd.^2)) - ...
        (t/T).*exp((-(x-OD).^2)./(2*sd.^2));
    k = k+1;
end
V = zeros(N,T);
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
end