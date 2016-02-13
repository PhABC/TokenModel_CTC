%% FigureCTC 
% script for all the figures

%% Stimuli
subplot(4,1,1); 
plot(sign(S.stimW)*S.stim(trial,:),'x','linewidth',1); hold on;
plot(zeros(S.T,1),'m');                      hold off

title('Stimuli');
ylabel('Token right - Token Left');
xlabel('time (ms)');
xlim([0,S.T])
ylim([-16,16])

%% Prefered direction VS non prefered
subplot(4,1,2)
plot(mean(M1(S.npref,:)),'r','linewidth',2); hold on
plot(mean(M1(S.pref, :)),'g','linewidth',2); hold off

set(gca,'Linewidth',3);
set(gca,'FontSize',40);
title('Prefered direction VS non prefered');
xlabel('time (ms)');
ylabel('activation');
xlim([0,S.T])
ylim([0,101])


%% Population 1 activity
subplot(4,1,3)
imagesc(M1)
title('Population 1 activity');
xlabel('time (ms)');
ylabel('Neurons');
xlim([0,S.T])
caxis([0,  S.beta ])

%% Population 2 activity
subplot(4,1,4)
imagesc(M2)

title('Population 2 activity');
xlabel('time (ms)');
ylabel('Neurons');
xlim([0,S.T])
caxis([0,  S.beta ])
hold off

drawnow;
