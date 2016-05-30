function [stim,idx,S] = ExtInputs(S,Stim)
%% External inputs
    %All the inputs not coming from the neural activity

%% Stimuli
%Unpacking fields
T=S.T; dt=S.dt; nbEx=S.nbEx; c=S.c;jumpT=S.jumpT; N = S.N;
Uslop=S.Uslop; onset=S.onset; Uw=S.Uw; Uori = S.Uori;
urgmax=S.urgmax;

%Matrix to expand stimuli wrt T with linspace
timeStimMat = zeros(15,T);
i = 1;

for bin = onset:jumpT:T
   if i == 1    
     timeStimMat(i,bin:bin+jumpT-1) = linspace(0,1,jumpT);
   elseif i<=15 && i>1
     timeStimMat(i,bin:bin+jumpT-1) = linspace(0,1,jumpT);
     timeStimMat(i-1,bin:bin+jumpT-1) = linspace(1,0,jumpT);
   else  
     timeStimMat(15,bin:end) = 1;
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
    		idx = randi(nstim,[nbEx,1]);  %Select nbEx random integer within nstim
	else    
    		lenC = cellfun('length',Stim.idx(S.c)); % Number of stim in each class
    		repD = ceil( ones(1,length(lenC))./lenC * max(lenC)); % Factor to reshape
    	
    		stUnpck  = Stim.idx(S.c); % Unpacking stimuli
            
    		% Reshaping for more similar number of trials
            for i = 1:length(stUnpck)
                stUnpck{i} = repmat(stUnpck{i},repD(i),1);
    		end
	  				
    		stimList = vertcat(stUnpck{:});             % List of stimuli
    		perms    = randi(length(stimList),nbEx,1);  % random stimuli
    		idx      = stimList(perms);                 % Selecting stimuli
	end
end
   

%Stimuli from c expanded wrt T
stimAll = Stim.stimRawR*timeStimMat;

%Choosing examples
stim = stimAll(idx,:);
%    
% %% Urgency   
% %Baseline urgency signal
%    Uend = (T-onset)*Uslop;       %Last point of linear function wrt T
%    urg  = zeros(nbEx,T);  
%    urg(:,1:onset)   = Uori;
%    urg(:,onset:end) = repmat(linspace(Uori,Uend,T-onset+1),... 
%                              [nbEx,1]);      % Urgency starts at stim onset
% 
% %Creating gauss distribution   
% if S.URtrial(1)
%    urg  = urg  + urg  .* repmat(S.URtrial(2)*randn(nbEx,1),1,T) + ...  % Random slopes
%           Uori + Uori .* repmat(S.URtrial(2)*randn(nbEx,1),1,T);    	             % Random origins
% end
% 
%    urg(urg<0) = 0;                                        % Del values < 0
% 
%    urg = urg*Uw;
%    urg(urg>urgmax) = urgmax; %max urg 
%     


%% Adding randomness to parameters
    
% Activation function
if S.RFsteep(1) ==1 
    S.STEEP = randParam(S.Fsteep,S.RFsteep(2),S.N);
else
    S.STEEP = repmat(S.Fsteep,N,1);
end

if S.RFshift(1) ==1
    S.SHIFT = randParam(S.Fshift,S.RFshift(2),S.N);
else
    S.SHIFT = repmat(S.Fshift,N,1);
end
