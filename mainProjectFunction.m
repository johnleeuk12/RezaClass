function mainProjectFunction(x0, A,B, pSteps, k1Steps, T, Thold, ...
    numAverage, Qx, Qy, H, C1, C2, L, Delta, goal, ind1, ind2)
%% computing Kalman gains and Optimal feedback gains


total_steps = k1Steps+pSteps;

for k = 1:total_steps
    G{k} = zeros(2,7);
end

for i =  1:10
    Kalman = calculateKalmanGains(x0,G,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B);
    G = calculateFeedbackGains(Kalman,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B,T,L);
%     pause
end

%% Intializing parameters for simulation
x_hat{1} = x0;
x{1} = x0;
phi_1 = randn(1,total_steps);
phi_2 = randn(1,total_steps);
% eta_1 = randn(1,pSteps);
% eta_2 = randn(1,pSteps);
eps_x = mvnrnd(zeros(size(x0)),Qx,total_steps).';
eps_y = mvnrnd(zeros(size(H*x0)),Qy,total_steps).';
y = zeros(2,total_steps);


%% Simulation
for k = 1: total_steps-1
    u{k} = -G{k}*x_hat{k};
    y(:,k) = H*x{k} + eps_y(:,k);
    x{k+1} = A*x{k} + B*(u{k} + C1*u{k}*phi_1(k) + C2*u{k}*phi_2(k)) + eps_x(:,k);
    x_hat{k+1} = A*x_hat{k} + A*Kalman{k}*(y(:,k)-H*x_hat{k}) + B*u{k};
    uu(:,k) = u{k}.';
end

figure
plot(y(1,:));
hold on
plot(y(2,:));
% plot(y(2,:)+y(1,:));
test = 1;
function Kalman = calculateKalmanGains(x0,G,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B)

% compute Kalman gain with given initial position x0 and Feedback Gain G

%seed and intializing parameters
S_x{1} = x0*x0.';
S_e{1} = zeros(7,7);
S_xe{1} = zeros(7,7);
total_steps = k1Steps+pSteps;

for k = 1:total_steps-1
    % computing the Kalman filter
    Kalman{k} = S_e{k}*H.'/(H*S_e{k}*H.'+Qy);
    S_e{k+1} = A*(eye(7)-Kalman{k}*H)*S_e{k}*A.' + Qx + B*C1*G{k}*S_x{k}*G{k}.'*C1.'*B.' ...
        + B*C2*G{k}*S_x{k}*G{k}.'*C1.'*B.';
    S_x{k+1} = (A-B*G{k})*S_x{k}*(A-B*G{k}).' + (A-B*G{k})*S_xe{k}*(A*Kalman{k}*H).' ...
        + A*Kalman{k}*H*S_xe{k}.'*(A-B*G{k}).' + A*Kalman{k}*H*S_e{k}*A.';
    S_xe{k+1} = (A-B*G{k})*S_xe{k}*(eye(7)-Kalman{k}*H).'*A.';
    
end

Kalman{total_steps} = S_e{total_steps}*H.'/(H*S_e{total_steps}*H.'+Qy);


function G = calculateFeedbackGains(Kalman,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B,T,L)
%compute the feedback gain with position x and the Kalman filter Kalman

total_steps = k1Steps+pSteps;
%seed 
for k = 1:total_steps
    W_x{k} = zeros(7,7);
    W_e{k} = zeros(7,7);
    w(k) = 0;
end
 
%Initialize at end time p = total_steps

W_x{total_steps} = H.'*T*H;

for k = total_steps-1: -1 : 1
    if k <=k1Steps
        T= zeros(2,2);
    end
    % computing feedback Gain G
    Cx = C1.'*B.'*W_x{k+1}*B*C1 + C2.'*B.'*W_x{k+1}*B*C2;
    Ce = C1.'*B.'*W_e{k+1}*B*C1 + C2.'*B.'*W_e{k+1}*B*C2;
    %D = 
    G{k} = (L + Cx + Ce + B.'*W_x{k+1}*B)\B.'*W_x{k+1}*A;
    
    W_x{k} = H.'*T*H + A.'*W_x{k+1}*A - A.'*W_x{k+1}*B*G{k};
    W_e{k} = A.'*W_x{k+1}*B*G{k} + (A-A*Kalman{k}*H).'*W_e{k+1}*(A-A*Kalman{k}*H);
    w(k) = trace(W_x{k+1}*Qx + W_e{k+1}*(Qx + A*Kalman{k}*Qy*Kalman{k}.'*A.'));
end

G{total_steps} = zeros(2,7);



















