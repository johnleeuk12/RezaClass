function hw23_jlee629()

%% Initializing Parameters

A = [1 0.001 0;-0.34 0.92 0.34; 0 0 0.9];
b = [0;0;0.1];
s = [1;0;0];
dt = 1e-3;

p = 120;
k = 50;
kappa = 0;


Z1 = zeros(p+k-1,p+k-1);
for n = 1:p+k-1
    Z1(n,n) = 2*kappa*s.'*(A^(p+k-n))*(b*b.')*(A^(p+k-n)).'*s;
end

Z2 = zeros(p+k-1,k);
for n = 1:p+k-1
    for i = 1:k
        Z2(n,i) = s.'*(A^(p+i-2-n))*b;
    end
end

Z3 = zeros(k,p+k-1);
for n = 1:p+k-1
    for i = 1:k
        Z3(i,n) = s.'*(A^(p+i-1-n))*b;
    end
end

Z4 = zeros(k,k);

Z = [Z1 Z2;Z3 Z4];

