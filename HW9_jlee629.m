function HW9_jlee629()

clear all

load('Kalman_HW8_dat.mat');

%% initializing Parameters
A{1} = [0.5 0; 0 0.3];
A{2} = eye(2);
A{3} = eye(2);
A{4} = eye(2);

epsX{1} = eye(2);
epsX{2} = [25 0; 0 25];
epsX{3} = eye(2);
epsX{4} = eye(2);

epsS{1} = eye(2);
epsS{2} = eye(2); 
epsS{3} = eye(2);
epsS{4} = [25 0; 0 25];

epsY{1} = 1.;
epsY{2} = 1.;
epsY{3} = 25;
epsY{4} = 1.;

N = length(y); 
H = [1 1];

%% initializing the system
X = zeros(2,N);
Xhat = zeros(2,N);
Yhat = zeros(1,N);
p = 1;
X(:,1) = A{1}*[0;1] + normrnd(zeros(2,1),sqrt([epsX{p}(1,1); epsX{p}(2,2)])); % this is X(n)
Xhat(:,1) = [0;1];
Yhat(1) = H*(Xhat(:,1) + normrnd(zeros(2,1),sqrt([epsX{p}(1,1); epsX{p}(2,2)]))) + normrnd(zeros(1,1), sqrt(epsY{p}));
P{1} = [1 0; 0 1];
% K = {};
K(:,1) = [0;0];
LMS(:,1) = [0;0];
%% computing the Kalman gain for different states. 
for n = 2:N
    if n<201
        p = 1;
    elseif n < 401
        p = 2;
    elseif n < 601
        p = 3;
    else
        p = 4;
    end
    X(:,n) = x(:,n);
%     LMS(:,n) = X(:,n)/(X(:,n).'*X(:,n)); %computing LMS
    K(:,n) = P{n-1}*H.'/(H*P{n-1}*H.' +H*(sqrt(X*X.').*epsS{p})*H.' + epsY{p}); %computing Kalman gain 
    Ptemp = (eye(2) - K(:,n)*H)*P{n-1};
    Xhat(:,n) =(eye(2)- K(:,n)*H)*Xhat(:,n-1) + K(:,n)*(H*X(:,n)+ H*normrnd(zeros(2,1),sqrt(abs(X(:,n)).*[epsX{p}(1,1); epsX{p}(2,2)])) + normrnd(zeros(1,1), sqrt(epsY{p})));
    P{n} = A{p}*Ptemp*A{p}.' + epsX{p};
    Yhat(n) = H*(Xhat(:,n));
end


%% figures

figure
plot(K(1,:))
title('Kalman Filter')
xlabel('iterations')
% hold on 
% plot(LMS(1,:))
figure
title('Estimate of X')
plot(Xhat(1,:))
hold on 
plot(Xhat(2,:))
xlabel('iterations')



figure
title('comparison between y and Yhat')
plot(Yhat)
hold on 
plot(y)
xlabel('iterations')
legend('Yhat','Y')


% 
% figure
% plot(y-Yhat)
% 
% figure
% 
% test = zeros(1,N);
% for n = 2:N
%     test(n) = sqrt((LMS(1,n)-K(1,n))^2 + (LMS(2,n)-K(2,n))^2);
% end
% plot(test)
% 
% 
% 
% 
