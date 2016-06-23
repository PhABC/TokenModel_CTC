function clusterPlot(X,norm,conds,comps)
% Will plot clusters of loading matrices of first 3 components
% combining PMd and M1 cells. 

if nargin == 1
    norm  = 1;
    conds = 1:6;
    comps = [1,2,4];
elseif nargin == 2
    conds = 1:6; 
    comps = [1,2,4];
elseif nargin == 3
    comps = [1,2,4];
end

if length(comps) ~= 3
    error('Need 3 components : comps should be a vector of length 3')
end

% X is the data
loads = pcaTok(X,norm,conds);
loads = loads.full;

% "Normalizing" loadings
if norm
    loadnorm = repmat(max(abs(loads)),size(loads,1),1);
    loads    = loads./loadnorm;
end 


c1 = comps(1);
c2 = comps(2);
c3 = comps(3);

%M1
% 1 & 2
figC1 = figure; hold on
subplot(2,3,1);
plot(loads(301:end,c1),loads(301:end,c2),'.r','MarkerSize',10)

xlim([-1,1])
ylim([-1,1])
xlabel([ 'Component ', num2str(c1) ])
ylabel([ 'Component ', num2str(c2) ])
drawnow

% 1 & 3
%figC2 = figure; hold on
subplot(2,3,2); hold on
plot(loads(301:end,c1),loads(301:end,c3),'.r','MarkerSize',10)

xlim([-1,1])
ylim([-1,1])
xlabel([ 'Component ', num2str(c1) ])
ylabel([ 'Component ', num2str(c3) ])
drawnow

% 2 & 3
% figC3 = figure; hold on
subplot(2,3,3); hold on
plot(loads(301:end,c2),loads(301:end,c3),'.r','MarkerSize',10)

xlim([-1,1])
ylim([-1,1])
xlabel([ 'Component ', num2str(c2) ])
ylabel([ 'Component ', num2str(c3) ])
drawnow

%% PMd

% 1 & 2
%figure(figC1); hold on
subplot(2,3,1); hold on
plot(loads(1:300,c1),loads(1:300,c2),'.b','MarkerSize',10)
drawnow

% 1 & 3
%figure(figC2); hold on
subplot(2,3,2); hold on
plot(loads(1:300,c1),loads(1:300,c3),'.b','MarkerSize',10)
drawnow

% 2 & 3
%figure(figC3); hold on
subplot(2,3,3); hold on
plot(loads(1:300,c2),loads(1:300,c3),'.b','MarkerSize',10)
drawnow




