function hw22_jlee629()

%% Parameter initialization

Pa = [9 49 100 145 190 230 270 340]; % deg
Pa = Pa/180*pi;
P = [cos(Pa); sin(Pa)];
lambda = 0.5;
m = 1.1;
Z = {};
f_ga = [1:2:360]/180*pi;
u_opt = zeros(length(f_ga),8);
f_g = [cos(f_ga);sin(f_ga)];
A = eye(8)*-1;
B = zeros(8,1);
%% optimization of parameter u
for n = 1:length(f_ga)
    u0 = rand(1,8);
    J = @(u)(f_g(:,n)-P*u.').'*(f_g(:,n)-P*u.') + lambda*sum(u.^m); % cost function
    u_opt(n,:) = fmincon(J,u0,A,B); %calculating optimal force u
    disp(n)
end

%% projecting u on f_g
for n = 1:length(f_ga)
    Z.x(n,:) = u_opt(n,:)*f_g(1,n);
    Z.y(n,:) = u_opt(n,:)*f_g(2,n);
end


figure
for t = 1:8
    scatter(Z.x(:,t),Z.y(:,t),10,'filled')
    hold on
end
plotv(P,'-')

