%% 

subplot(2,1,1)
%plot(zeros(S.N,1))
plot(fct([0:1:100],S.Fsteep(1),S.Fshift(1)))

subplot(2,1,2)
plot(fct([0:1:100],S.Fsteep(2),S.Fshift(2)))

