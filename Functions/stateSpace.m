function stateSpace(X,comps,conds)
% StateSpace will plot the neural activity for each cond projected on the 
% first 3 components in comps found by PCA. 

if nargin == 1
    comps = 1:6;
    conds = 1:6;
end

% Stack cells wrt fields
Xnet = cell(length(X),1);
for net = 1:length(X)
    Xnet{net} = stackStruct(X{net});  % stackStruct is adapted to this FR format
end
Xfr = stackStruct(Xnet);
    
% Symetrise cells for left direction
% X2fr = Xfr;
X2fr = symCells(Xfr);
X2fr = structfun( @(X) X(:,:,conds), X2fr, ...
                   'UniformOutput', false);

N = size(X2fr.PMd,1);
T = size(X2fr.PMd,2);
C = size(X2fr.PMd,3);

% Normalizing for max FR
X2fr = structfun( @(X) sqrt(X),X2fr, ... 
                  'UniformOutput', false);

fn = fieldnames(X2fr); 
for f = 1:length(fn)
    Xresh.(fn{f}) = reshape(X2fr.(fn{f})(:,5:T-7,:),N,(T-7-4)*C);
    [loads.(fn{f}),scores.(fn{f}),latent.(fn{f})] = pca(Xresh.(fn{f})'); 
end
 
%% Building components


for f = 1:length(fn)
    for comp = comps
        a = 0;
        for cond = conds
            a = a+1;
            build.(fn{f})(:,comp,a) =  sum( diag(loads.(fn{f})(:,comp))*...
                                            X2fr.(fn{f})(:,5:T-7,cond) )';
        end
    end
end


%% Plotting

color = [ 'b', 'b', 'r', 'r', 'g', 'g','y','y' ];
% color = 'k';

% Neural state space
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

% First 6 PCs
for f = 1:length(fn)
    figure; hold on;    
    for c = 1:6
        for cond = 1:a
            subplot(2,3,c); hold on
            plot(build.(fn{f})(:,c,cond),'color', color(cond))
            
            title( [fn{f},' Component ', num2str(c), ' ~ R² = '...
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

fn = fieldnames(X); 

N  = size(X.(fn{1}),1);
% T  = size(X.(fn{1}),2);
C  = size(X.(fn{1}),3);
    

for f = 1:length(fn)-3 
    
    for c = 1:C*2
        if mod(c,2) 
           X2.(fn{f})(:,:,c) = X.(fn{f})(:,:,(c-1)/2+1);
        else
           X2.(fn{f})(:,:,c) = circshift(X.(fn{f})(:,:,(c-2)/2+1),N/2);
           %X2.(fn{f})(:,:,c)  = X.(fn{f})(:,:,(c-2)/2+1);
        end
    end
     
end

