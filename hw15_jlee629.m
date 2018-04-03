function hw15_jlee629()
clear all


Data = load('classification_IRLS.dat');
N = length(Data); 
X = [ones(N,1) Data(:,2) Data(:,3)];

w = [0 0 0];

Q = zeros(N,N);
for it = 1:10
    for n = 1:N
        q(n) = 1/(1+exp(-w(it,:)*X(n,:).'));
        Q(n,n) = q(n)*(1-q(n));
    end
    w(it+1,:) = w(it,:)+ (inv(X.'*Q*X)*X.'*(Data(:,1)-q.')).';
    err(it) = immse(Data(:,1),q.');
end

figure
for cls = 0:1
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
end

syms x1 x2
eqn = w(end,:)*[1 x1 x2].';

fimplicit(eqn)

axis([-7,3,-1,7]);

hold off

figure
plot(err) 