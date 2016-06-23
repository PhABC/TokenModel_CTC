%% FigureCTC 
% script for all the figures

figure(99)

time = linspace(0,S.T*S.dt,S.T);

%% Stimuli
subplot(4,1,1); 
plot(time,sign(S.stimW)*stimTrial(trial,:),'x','linewidth',1); hold on;
plot(zeros(S.T,1),'m');                      
plot(S.onset*S.dt,-16:001:16,'xk'); hold off

title('Stimuli');
ylabel('Token right - Token Left');
xlabel('time (ms)');
xlim([0,S.T*S.dt])
ylim([-16,16])

%% Prefered direction VS non prefered
subplot(4,1,2)
if ~commit(trial)
    plot(time,mean(PMd_(npref,:)),'r','linewidth',2); hold on
    plot(time,mean(PMd_(pref ,:)),'g','linewidth',2);
else
    timeTrial = floor(1:min(abs(commit(trial)/S.dt)+S.aftcmt,S.T));
    plot(time(timeTrial), mean(PMd_(npref,timeTrial)),'r','linewidth',2); hold on
    plot(time(timeTrial), mean(PMd_(pref, timeTrial)),'g','linewidth',2);
end
plot(S.onset*S.dt,0:01:S.beta,'xk'); 
plot(abs(commit(trial)/S.dt),0:2:S.beta,'xb'); hold off

% set(gca,'Linewidth',3);
% set(gca,'FontSize',12);
title('Prefered direction VS non prefered');
xlabel('time (ms)');
ylabel('activation');
xlim([0,S.T*S.dt])
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
