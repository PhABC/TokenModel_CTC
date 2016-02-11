function [w1,w2] = WeightMatrix(N)
%% Weight matrix creation
    
w = zeros(N);

for i = 1:N
    w(i,max(i-3,1):i) = linspace(0,0.4,length(max(i-3,1):i));
    w(i,i:min(i+3,N)) = linspace(0.4,0,length(i:min(i+3,N)));
end

w1 = w+randn(N).*max(max(w)).*0.01;
% w2 = w+randn(N).*max(max(w)).*0.01; %recurrent connections
w2 = zeros(N); %FF connections
 

asdffd 
end

