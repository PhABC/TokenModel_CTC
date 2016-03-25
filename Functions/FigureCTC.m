%% FigureCTC 
% script for all the figures

%% Stimuli
subplot(4,1,1); 
plot(sign(S.stimW)*stimTrial(trial,:),'x','linewidth',1); hold on;
plot(zeros(S.T,1),'m');                      
plot(S.onset,-16:001:16,'xk'); hold off


title('Stimuli');
ylabel('Token right - Token Left');
xlabel('time (ms)');
xlim([0,S.T])
ylim([-16,16])

%% Prefered direction VS non prefered
subplot(4,1,2)
if ~commit(trial)
    plot(mean(PMd_(npref,:)),'r','linewidth',2); hold on
    plot(mean(PMd_(pref ,:)),'g','linewidth',2);
else
    plot(mean(PMd_(npref,1:min(abs(commit(trial))+S.aftcmt/S.dt,S.T))),'r','linewidth',2); hold on
    plot(mean(PMd_(pref, 1:min(abs(commit(trial))+S.aftcmt/S.dt,S.T))),'g','linewidth',2);
end
plot(S.onset,0:001:S.beta,'xk'); hold off

set(gca,'Linewidth',3);
set(gca,'FontSize',40);
title('Prefered direction VS non prefered');
xlabel('time (ms)');
ylabel('activation');
xlim([0,S.T])
ylim([0,101])


%% Population PMd activity
subplot(4,1,3)
imagesc(circshift(PMd_,ceil(S.N/4),1)); hold on
plot(S.onset,0:001:S.N,'xk'); hold off

title('PMd population  activity');
xlabel('time (ms)');
ylabel('Neurons');
xlim([0,S.T])
caxis([0,  S.beta ])

%% Population M1 activity
subplot(4,1,4)
imagesc(circshift(M1_,ceil(S.N/4),1)); hold on
plot(S.onset,0:001:S.N,'xk'); hold off

title('M1 population activity');
xlabel('time (ms)');
ylabel('Neurons');
xlim([0,S.T])
caxis([0,  S.beta ])
hold off

drawnow;
