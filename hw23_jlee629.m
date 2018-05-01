function hw23_jlee629()

%% Initializing Parameters

A = [1 0.001 0;-0.34 0.92 0.34; 0 0 0.9];
b = [0;0;0.1];
s = [1;0;0];
% dt = 1e-3;

p = 120;
k = 50;
kappa = 0.0005;

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

c=1e-13;
Z = Z+c*eye(length(Z));

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
X = zeros(3,p);
noise = randn(length(U),1)*kappa^2;

for t = 1:p-1
    X(:,t+1) = A*X(:,t)+b*(U(t)+noise(t)*U(t)^2);
end


subplot(1,2,1)
plot(X(1,:))
hold on
subplot(1,2,2)
plot(X(2,:))
hold on

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


