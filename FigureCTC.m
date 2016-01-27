%% FigureCTC 
% script for all the figures

figure;
hold on;
plot(mean(M1(npref,:)),'r','linewidth',4);
plot(mean(M1(pref,:)),'g','linewidth',4);
set(gca,'Linewidth',3);
set(gca,'FontSize',40);
xlabel('time (ms)');
ylabel('activation');