function mainProjectFunction(x0, A,B, pSteps, k1Steps, T, Thold, ...
    numAverage, Qx, Qy, H, C1, C2, L, Delta, goal, ind1, ind2)
%% Initializing certain parameters

for k = 1:pSteps
    G{k} = zeros(2,7);
end

Ts = {};
for k = 1:pSteps
    if k < k1Steps
        Ts{k} = zeros(2,2);
    else
        Ts{k} = T;
    end
end

%% computing Kalman gains and Optimal feedback gains
for i =  1:10
    Kalman = calculateKalmanGains(x0,G,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B);
    G = calculateFeedbackGains(Kalman,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B,Ts,L);
end

%% Simulation
for n = 1:numAverage
    % Intializing parameters for simulation
    x_hat{1} = x0;
    x{1} = x0;  
    y = zeros(2,pSteps);
    
    phi_1 = randn(1,pSteps);
    phi_2 = randn(1,pSteps);
    eps_x = mvnrnd(zeros(length(Qx),1),Qx,pSteps).';
    eps_y = mvnrnd(zeros(length(Qy),1),Qy,pSteps).';
    eye_posi(n,1) = 0;
    head_posi(n,1) = 0;
    
    %simulation
    for k = 1: pSteps-1
        u{k} = -G{k}*x_hat{k};
        y(:,k) = H*x_hat{k} + eps_y(:,k) ; %
        x{k+1} = A*x{k} + B*u{k} + B*(C1*u{k}*phi_1(k) + C2*u{k}*phi_2(k)) + eps_x(:,k) ; %
        x_hat{k+1} = A*x_hat{k} + A*Kalman{k}*(y(:,k)-H*x_hat{k}) + B*u{k};
        uu(:,k) = u{k}.';
        eye_posi(n,k+1) = x{k+1}(1)/pi*180;
        head_posi(n,k+1) = x{k+1}(4)/pi*180;
        %hold time
        if Thold ~= 0
            if k < Thold
                x{k+1}(4) = 0;
                x{k+1}(5) = 0;
                x_hat{k+1}(5) = 0;
                x_hat{k+1}(4) = 0;
            end
        end
    end
end


eye_mean = mean(eye_posi,1);
head_mean = mean(head_posi,1);
timebin = Delta:Delta:pSteps*Delta;
%% plotting figures
figure
plot(timebin,eye_mean,'LineWidth',2,'DisplayName','eye');
hold on
plot(timebin,head_mean,'LineWidth',2,'DisplayName','head');
plot(timebin,eye_mean+head_mean,'LineWidth',2,'DisplayName','gaze');
plot(timebin,40*ones(1,length(timebin)),'--k','LineWidth',2)

legend('show')
xlabel('time (s)')
ylabel('position (deg)')
title(['Figure n. ' num2str(ind2) ' simulation average']) 
grid('on')


hold off

figure
plot(timebin,eye_posi(1,:),'LineWidth',2,'DisplayName','eye');
hold on
plot(timebin,head_posi(1,:),'LineWidth',2,'DisplayName','head');
plot(timebin,eye_posi(1,:)+head_posi(1,:),'LineWidth',2,'DisplayName','gaze');
plot(timebin,40*ones(1,length(timebin)),'--k','LineWidth',2)

legend('show')
xlabel('time (s)')
ylabel('position (deg)')
title(['Figure n. ' num2str(ind1) ' single simulation']) 
grid('on')
hold off
test = 1;
function Kalman = calculateKalmanGains(x0,G,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B)

% compute Kalman gain with given initial position x0 and Feedback Gain G

%seed and intializing parameters
S_x{1} = x0*x0';
S_e{1} = zeros(7,7);
S_xe{1} = zeros(7,7);


for k = 1:pSteps-1
    % computing the Kalman filter
    Kalman{k} = S_e{k}*H' * inv((H*S_e{k}*H'+Qy));
    S_e{k+1} = A*(eye(7)-Kalman{k}*H)*S_e{k}*A' + Qx + B*C1*G{k}*S_x{k}*G{k}'*C1'*B' ...
        + B*C2*G{k}*S_x{k}*G{k}.'*C1.'*B.';
    S_x{k+1} = (A-B*G{k})*S_x{k}*(A-B*G{k})' + (A-B*G{k})*S_xe{k}*(A*Kalman{k}*H)' ...
        + A*Kalman{k}*H*S_xe{k}'*(A-B*G{k})' + A*Kalman{k}*H*S_e{k}*A';
    S_xe{k+1} = (A-B*G{k})*S_xe{k}*(eye(7)-Kalman{k}*H)'*A';
    
end

Kalman{pSteps} = S_e{pSteps}*H.'/(H*S_e{pSteps}*H.'+Qy);


function G = calculateFeedbackGains(Kalman,pSteps, k1Steps,Qx, Qy, H, C1, C2,A,B,Ts,L)
%compute the feedback gain with position x and the Kalman filter Kalman


%seed 
for k = 1:pSteps
    W_x{k} = zeros(7,7);
    W_e{k} = zeros(7,7);
%     w(k) = 0;
end
 
%Initialize at end time p = total_steps

W_x{pSteps} = H.'*Ts{pSteps}*H;


for k = pSteps-1: -1 : 1
    % computing feedback Gain G
    Cx = C1.'*B.'*W_x{k+1}*B*C1 + C2.'*B.'*W_x{k+1}*B*C2;
    Ce = C1.'*B.'*W_e{k+1}*B*C1 + C2.'*B.'*W_e{k+1}*B*C2;
    De = zeros(size(W_x{k+1}));
    G{k} = (L + Cx + Ce + B' * W_x{k+1} * B) \ (B' * W_x{k+1} * A);
    
    W_x{k} = H' * Ts{k} * H + A' * W_x{k+1} * A - G{k}' * B' * W_x{k+1} * A + De;
    W_e{k} = G{k}.'*B.'*W_x{k+1}*A + (A-A*Kalman{k}*H).'*W_e{k+1}*(A - A * Kalman{k} * H);
end

G{pSteps} = zeros(2,7);



















