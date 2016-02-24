function [PMD,M1,commit]  = trialLoop(S,net)

    commit = zeros(1,S.nbEx);
    PMD    = cell(1,S.nbEx);
    M1     = cell(1,S.nbEx);
    for trial = 1:S.nbEx

        % stimuli preferences are also assigned in the following fct
        [ St, V, U, pref, npref ] = TrialInput(S,net,trial);
        [ PMD_,M1_,commit(trial)] = CTCsim(S,net,St,V,U,pref,npref); 

        if S.printDec
            
            %Number of tokens (15 is the max)
            nbTokens =  floor((abs(commit(trial))-S.onset)/50)+1;
            if nbTokens > 15; nbTokens = 15; end
        
            fprintf('      Trial :   %d of %d    ~    Decision :   ',trial,S.nbEx)
            
            if commit(trial) > 0
                fprintf('Right (tokens: %d)\n', nbTokens);
            elseif commit(trial) < 0
                fprintf('Left  (tokens: %d)\n', nbTokens);
            else 
                fprintf('--- \n')
            end
            
        end

        %Figures
        if S.plotting 
           FigureCTC
        end
        
        %Only outputting values of prefered neurons
        PMD{trial} = PMD_(pref,1:min(abs(commit(trial))+S.aftcmt,S.T));
        M1{trial}  =  M1_(pref,1:min(abs(commit(trial))+S.aftcmt,S.T));      
    end

