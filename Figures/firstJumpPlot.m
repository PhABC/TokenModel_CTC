% Will plot single cell activity

X{1}{1} = firstJumpFormat(100,FR{1}{1},Stim,idxStim);
X{2}{1} = firstJumpFormat(100,FR{2}{1},Stim,idxStim);

% Concatenating

ncells = length(X);                        

C = X{1}{1};

fn = fieldnames(C);                                                                                                        
   
for c = 2:ncells                                                                                                                       
    for f = 1:length(fn)                                                                                                         
        C.(fn{f}) = vertcat(C.(fn{f}),X{c}{1}.(fn{f}));                                    
    end                                           
end 

pref  = zeros(212,1);
pref(1:26)    = 1;
pref(80:132)  = 1; 
pref(186:end) = 1;
pref = logical(pref);
npref = ~pref;

%% PMd

figure
%pref
plot(mean(C.PMd(pref,:,1)),'b','LineWidth',1.3); hold on;
plot(mean(C.PMd(pref,:,2)),'r','LineWidth',1.3)
plot(mean(C.PMd(pref,:,3)),'g','LineWidth',1.3)
%npref
plot(mean(C.PMd(npref,:,1)),'.b','MarkerSize',5)
plot(mean(C.PMd(npref,:,2)),'.r','MarkerSize',5)
plot(mean(C.PMd(npref,:,3)),'.g','MarkerSize',5)
ylim([5,100])
title('PMd')

%% M1

figure
%pref
plot(mean(C.M1(pref,:,1)),'b','LineWidth',1.3); hold on;
plot(mean(C.M1(pref,:,2)),'r','LineWidth',1.3)
plot(mean(C.M1(pref,:,3)),'g','LineWidth',1.3)
%npref
plot(mean(C.M1(npref,:,1)),'.b','MarkerSize',5)
plot(mean(C.M1(npref,:,2)),'.r','MarkerSize',5)
plot(mean(C.M1(npref,:,3)),'.g','MarkerSize',5)
ylim([5,100])
title('M1')

