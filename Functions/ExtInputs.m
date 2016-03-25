function [ urg,stim,SG,idx ] = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli
%Unpacking fields
T=S.T; dt=S.dt; nbEx=S.nbEx; c=S.c;jumpT=S.jumpT; 
Uslop=S.Uslop; onset=S.onset; Uw=S.Uw; Uori = S.Uori;
urgmax=S.urgmax; exRan=S.exRan; exNum=S.exNum;

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

%Stimuli index e
if ~S.exRan
	idx = ones(S.nbEx,1)*S.exNum;
else
	if ~c
    		nstim   = Stim.nbStimR;       % Number of stim of chosen 
    		idx = randi(nstim,[nbEx,1]); %Select nbEx random integer within nstim
	else    
    		lenC = cellfun('length',Stim.idx(S.c));
    		repD = ceil( ones(1,length(lenC))./lenC * max(lenC));    
    	
    		stUnpck  = Stim.idx(S.c);
    		for i = 1:length(stUnpck)
			stUnpck{i} = repmat(stUnpck{i},repD(i),1);
    		end
	  				
    		stimList = vertcat(stUnpck{:}); 
    		perms    = randi(length(stimList),nbEx,1);
    		idx      = stimList(perms);
	end
end

   
stimAll = Stim.stimRawR*timeLogiMat;   %Stimuli from c expanded wrt T
   
   
%Choosing examples

stim = stimAll(idx,:);
   
%% Urgency   
%Baseline urgency signal
   Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
   urg  = zeros(nbEx,T);  
   urg(:,onset:end) = repmat(linspace(0,Uend,T-onset+1),... 
                             [nbEx,1]);      % Urgency starts at stim onset
%    urg  = urg/max(urg(:,end));                             % Normalized


%Creating gauss distribution   
if S.Urand
   urg  = urg  + urg .* abs(repmat(randn(nbEx,1)*.1,1,T)) + ...  % Random slopes
          Uori + repmat(0.1*randn(nbEx,1),1,T);    	         % Random origins
end
%    urg  = urg/max(urg(:,end));                            % Normalized
   urg(urg<0) = 0;                                        % Del values < 0

   urg = urg*Uw;
   urg(urg>urgmax) = urgmax; %max urg 

%% Slow noise
    tauSG = 500;                          % time constant of slow noise in ms
    noise = randn(nbEx,T/(tauSG/dt)+1);
    
    step = tauSG/dt; 
    
    SG = zeros(nbEx,T);
    i = 1; 
    
    %Expanding slow noise
    for t = 1:step:T
        SG(:,t:t+step-1) = linspaceMat(noise(:,i),noise(:,i+1),step);
        i=i+1;
    end
    
    %Smoothing out noise
    for ex = 1:nbEx
        SG(ex,:) = smooth(SG(ex,:),100/dt);
    end
    
end
