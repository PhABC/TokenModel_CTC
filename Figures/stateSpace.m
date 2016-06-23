function [loads,build,scores,latent] = stateSpace(X,norm,conds,comps)
% StateSpace will plot the neural activity for each cond projected on the 
% first 3 components in comps found by PCA. 

plot3D = 0; % 0 : will not plot 3d neural state space | 1 : will plot

if nargin == 1
    norm  = 0;
    comps = 1:6;
    conds = 1:6;
elseif nargin == 2
    comps = 1:6;
    conds = 1:6;
elseif nargin == 3
    comps = 1:6;
end

[loads,build,scores,latent] = pcaTok(X,norm,conds,comps);
fn = fieldnames(loads);


%% Plotting
a = length(conds);
color = [ 'b', 'b', 'r', 'r', 'g', 'g','y','y' ];
% color = 'k';

% Neural state space
if plot3D
	for f = 1:length(fn)
    	figure; hold on;
    	title( [fn{f},' Neural State Space'])
    	for cond = 1:a
        	plot3(build.(fn{f})(:,3,cond),...
              	build.(fn{f})(:,2,cond),...
              	build.(fn{f})(:,1,cond),...
              	'color', color(cond)); hold on
    	end
    	xlabel(['PC3 ~ ', num2str((latent.(fn{f})(3)/...
             	sum(latent.(fn{f})))*100), '% variance explained']);
    	ylabel(['PC2 ~ ', num2str((latent.(fn{f})(2)/...
             	sum(latent.(fn{f})))*100),'% variance explained']);
    	zlabel(['PC1 ~ ', num2str((latent.(fn{f})(1)/...
             	sum(latent.(fn{f})))*100),'% variance explained']);
	end
end

% First 6 PCs
for f = 1:length(fn)
    figure('Name', fn{f}); hold on;    
    for c = 1:6
        for cond = 1:a
            if  ~mod(cond,2) 
		mark = '--';
	    else
		mark = '-';
	    end
	    subplot(2,3,c); hold on
            plot(build.(fn{f})(:,c,cond),mark,'color', color(conds(cond)),'LineWidth',1.5)
            
            title([' Dim ', num2str(c), '    ' ...
                    num2str((latent.(fn{f})(c)/sum(latent.(fn{f})))*100)])
            xlabel('Time')
        end
    end
end



function C = stackStruct(X)
% Will concatenante same fields in diffenret cells
% Cells in X must contain the same fields and dimensions

ncells = length(X);

if ~isstruct(X)
    C = X{1};
else
    C = X;
end
fn = fieldnames(C); 

for c = 2:ncells
    for f = 1:length(fn)-3
        C.(fn{f}) = vertcat(C.(fn{f}),X{c}.(fn{f}));
    end
end


function X2 = symCells(X)
% Will symmetrises cell by shifting cells to replicated opposite stimuli.

rmfields = {'Idx','Com','nbTokens'};
X        = rmfield(X,rmfields);
 
fn       = fieldnames(X); 
nfields  = length(fn);

N  = size(X.(fn{1}),1);
% T  = size(X.(fn{1}),2);
C  = size(X.(fn{1}),3);

for f = 1:nfields 
    
    for c = 1:C*2
        if mod(c,2) 
           X2.(fn{f})(:,:,c) = X.(fn{f})(:,:,(c-1)/2+1);
        else
           X2.(fn{f})(:,:,c) = circshift(X.(fn{f})(:,:,(c-2)/2+1),N/2);
           %X2.(fn{f})(:,:,c)  = X.(fn{f})(:,:,(c-2)/2+1);
        end
    end
     
end


