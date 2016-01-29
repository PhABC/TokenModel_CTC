function V = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli

%Unpacking fields
   T=S.T;  nbEx=S.nbEx; c=S.c;jumpT=S.jumpT; bias=S.bias; Uslop=S.Uslop;
   onset=S.onset; stimW=S.stimW; Utype=S.Utype; Uw=S.Uw;  Uori = S.Uori;
     
%Initializing
   nstim   = Stim.nbStim{c};               %Number of stim of chosen c
   stimRaw = Stim.stimRawR(Stim.idx{c},:); %Trials foc condition c
   
%Logic matrix to expand stimuli wrt T
   timeLogiMat = zeros(15,T);
   i = 1;
   
   for bin = onset:jumpT:T
       if i < 15
         timeLogiMat(i,bin:bin+jumpT-1) = 1;
       else
         timeLogiMat(i,bin:end) = 1;
         break;
       end
       i = i+1;
   end
   
   stimAll = stimRaw*timeLogiMat;   %Stimuli from c expanded wrt T
   
%Choosing examples
   idx     = randi(nstim,[nbEx,1]); %Select nbEx random integer within nstim
   stim    = stimAll(idx,:);
   
  
%% Urgency   

%Baseline urgency signal
   Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
   urg  = zeros(nbEx,T);  
   urg(:,onset:end) = repmat(linspace(0,Uend,T-onset+1),... 
                             [nbEx,1]); %Urgency starts at stim onset
   urg  = urg/max(urg(:,end));                             % Normalized
   
%Creating gauss distribution   
   urg  = urg .* abs(repmat(randn(nbEx,1)+Uslop,1,T))...  % Random slopes
               + repmat(0.1*randn(nbEx,1)+Uori,1,T);      % Random origins
   urg  = urg/max(urg(:,end));                            % Normalized
   urg(urg<0) = 0;                                        % Del values < 0

   
   if Utype == 1
       %Additive urgency
        V = stim*stimW  + urg*Uw + bias;

   elseif Utype ==2
       %Multiplicative urgency ~~~ NOT MULTIPLIYING ALL INPUTS CURRENTLY
        V = stim*stimW .* urg*Uw + bias;
       
   end
 
    
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