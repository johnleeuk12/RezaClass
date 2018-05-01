function hw23_jlee629()

%% Initializing Parameters

A = [1 0.001 0;-0.34 0.92 0.34; 0 0 0.9];
b = [0;0;0.1];
s = [1;0;0];
% dt = 1e-3;

p = 120;
k = 50;
kappa = 0.00006;

g = pi/5;
x0 = [0;0;0];

Z1 = zeros(p+k-1,p+k-1);
for n = 1:p+k-1
    Z1(n,n) = 2*(kappa^2)*s.'*(A^(p-n))*(b*b.')*(A^(p-n)).'*s;
end

Z2 = zeros(p+k-1,k);
for n = 1:p+k-1
    for i = 1:k
        Z2(n,i) = s.'*(A^(p+i-1-n))*b;
        if (p+i-1-n)<0
            Z2(n,i) = 0;
        end
    end
end

Z3 = zeros(k,p+k-1);
for n = 1:p+k-1
    for i = 1:k
        Z3(i,n) = s.'*(A^(p+i-1-n))*b;
        if (p+i-1-n)<0
            Z3(i,n) = 0;
        end
    end
end

Z4 = zeros(k,k);

Z = [Z1 Z2;Z3 Z4];

%% computing optimal motor commmand U
c=1e-13; %regularization factor.
Z = Z+c*eye(length(Z)); %Regularizing Z to make it inversible

Y = zeros(1,p+k-1);
Y1 = zeros(1,k);

for i = 1:k
    Y1(1,i) = g-s.'*(A^(p+i))*x0; 
end

% for z1 = 1:219    
%     for z2 = 1:219
%         if Z(z1,z2)==0
%             Z(z1,z2) = 0.1;
%         end
%     end
% end
Y = [Y Y1].';

U = inv(Z)*Y;
U = [0; U];
X = zeros(3,p+k);
%% Computing the trajectory X 
noise = randn(length(U),1)*kappa^2;

for t = 1:p+k-1
    X(:,t+1) = A*X(:,t)+b*(U(t)+noise(t)*U(t)^2);
end

subplot(2,2,1)
plot(U(1:p+k),'Linewidth',2);
hold on
subplot(2,2,3)
plot(X(1,:),'Linewidth',2,'DisplayName','\kappa')
title('motor commands')
xlabel('steps')
hold on
subplot(2,2,3)
plot([1 p+k],[g g]);
hold on
title('Position')
xlabel('steps')
subplot(2,2,4)
plot(X(2,:),'Linewidth',2)
title('Velocity ')
xlabel('steps')
hold on

    
