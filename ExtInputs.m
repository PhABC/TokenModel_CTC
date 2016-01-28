function V = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli

%Unpacking fields
   N=S.N; T=S.T;  nbEx=S.nbEx; c=S.c; hnorm=S.hnorm; jumpT=S.jumpT;
   start=S.start; stimW=S.stimW;
   
%Initializing
   nstim   = Stim.nbStim{c};              %Number of stim of chosen c
   stimRaw = Stim.stimRawR(Stim.idx{c},:); %Trials foc condition c
   
%Logic matrix to expand stimuli wrt T
   timeLogiMat = zeros(15,T);
   i = 1;
   
   for bin = start:jumpT:T
       if i < 15
         timeLogiMat(i,bin:bin+jumpT-1) = 1;
       else
         timeLogiMat(i,bin:end) = 1;
         break;
       end
       i = i+1;
   end
   
   stimAll = stimRaw*timeLogiMat; %Stimuli from c expanded wrt T
   
%Choosing examples
   idx   = randi(nstim,[nbEx,1]); % Select nbEx random integer within nstim
   stim  = stimAll(idx,:)*stimW;
   
   %Creating stim for simulation
   
  asdfasdf
       
  
   
   
   
    
%% Simple stimuli JPT    
% %easy, non-misleading stimulus
% stim = zeros(360,T);
% stim(PD,1:T) = 1:T;
% x = 1:360;
% k=1;
% for t = 1:T
%     stim(:,t) =  (t/T).*exp((-(x-PD).^2)./(2*sd.^2)) - ...
%         (t/T).*exp((-(x-OD).^2)./(2*sd.^2));
%     k = k+1;
% end
%     
% V = zeros(N,T);
% V = hnorm' * stim;
% V = V+repmat([1:T]./1000,N,1);  %urgency
% % imagesc(V)
% % return;
% 
% %misleading stimulus
% % PD = 270;
% % OD = 90;
% % stim = zeros(360,T);
% % stim(PD,1:T) = 1:T;
% % x = 1:360;
% % for t = 1:T-100
% %     stim(:,t) =  (t/T).*exp((-(x-OD).^2)./(2*sd.^2)) - ...
% %         (t/T).*exp((-(x-PD).^2)./(2*sd.^2));
% %     k = k+1;
% % end
% % 
% % for t = T-100:T
% %     stim(:,t) =  (t/T).*exp((-(x-PD).^2)./(2*sd.^2)) - ...
% %         (t/T).*exp((-(x-OD).^2)./(2*sd.^2));
% %     k = k+1;
% % end
% % V = hnorm' *stim;
% % imagesc(V)
% % return;
end