%% Parameter specification
clc; clear; close all;

% Time constants
tau1eye = 0.224;    % seconds
tau2eye = 0.013;    % seconds
alpha1eye = 0.004;  % seconds
tau1head = 0.270;   % seconds
tau2head = 0.15;    % seconds
alpha1head = 0.010; % seconds

% Covariances matrices for state and motor updates
Qx = 10^(-5) * eye(7);
Qy = 0.01 * eye(2);

% The goal direction
goal = 40 * pi/ 180;

% The H matrix
H = [-1,0,0,-1,0,0,1; 1,0,0,0,0,0,0];

% Specify the length of each time step
Delta = 0.002;  % seconds

% Specify the total simulation period
p = 0.3; % seconds

% The time at which we want our gaze to are at the target
k1 = 0.1; % seconds

% The number of discrete steps to reach time point k1
k1Steps = k1 / Delta;

% The number of time steps in the simulation
pSteps = p / Delta;

% The initial state vector
x0 = [0;0;0;0;0;0;goal]; 

% Specify the value of the T matrix. Note that T has this value only after
% k1 time has elapsed. It is zero before.
T =  [10^5,0; 0,300];

% Specify the matrix that determines the motor commands cost
L = [10,0;0,10];

% Number of simulations to average
numAverage = 50;

kalmanGains = cell(pSteps,5);
feedbackGains = cell(pSteps,5);

%% Generate A and B matrices
[A, B] = makeAandB(tau1eye, tau2eye, alpha1eye, tau1head, tau2head, alpha1head , Delta);

%% Case 1: No Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [0,0; 0,0];
C2 = [0,0; 0,0];

% The length of the holding period
Thold = 0;
figNum = 1;
figName = 'No Signal Dependent Noise, No Hold';
% Perform the simulations and plot the results
[Ks, Gs, eye_mean, head_mean, gaze_mean] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, figNum, figName);

kalmanGains(:,1) = reshape(mat2cell(Ks, size(Ks,1), size(Ks,2), ones(1,size(Ks,3))), 1, [])';
feedbackGains(:,1) = reshape(mat2cell(Gs, size(Gs,1), size(Gs,2), ones(1,size(Gs,3))), 1, [])';

%% Case 2: Small Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0   ,0; 0,0.01];

% The length of the holding period
Thold = 0;
figNum = 2;
figName = 'Small Signal Dependent Noise, No Hold';
% Perform the simulations and plot the results
[Ks, Gs, eye_mean, head_mean, gaze_mean] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, figNum, figName);

kalmanGains(:,2) = reshape(mat2cell(Ks, size(Ks,1), size(Ks,2), ones(1,size(Ks,3))), 1, [])';
feedbackGains(:,2) = reshape(mat2cell(Gs, size(Gs,1), size(Gs,2), ones(1,size(Gs,3))), 1, [])';

%% Case 3: Large Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [2.0 ,0; 0,0];
C2 = [0   ,0; 0,2.00];

% The length of the holding period
Thold = 0;
figNum = 3;
figName = 'Large Signal Dependent Noise, No Hold';
% Perform the simulations and plot the results
[Ks, Gs, eye_mean, head_mean, gaze_mean] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, figNum, figName);

kalmanGains(:,3) = reshape(mat2cell(Ks, size(Ks,1), size(Ks,2), ones(1,size(Ks,3))), 1, [])';
feedbackGains(:,3) = reshape(mat2cell(Gs, size(Gs,1), size(Gs,2), ones(1,size(Gs,3))), 1, [])';

%% Case 4: Small Signal Dependent Noise, Short Head Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0   ,0; 0,0.01];

% The length of the holding period
Thold = 0.050 / Delta;
figNum = 4;
figName = 'Small Signal Dependent Noise, Short Head Hold';
% Perform the simulations and plot the results
[Ks, Gs, eye_mean, head_mean, gaze_mean] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, figNum, figName);

kalmanGains(:,4) = reshape(mat2cell(Ks, size(Ks,1), size(Ks,2), ones(1,size(Ks,3))), 1, [])';
feedbackGains(:,4) = reshape(mat2cell(Gs, size(Gs,1), size(Gs,2), ones(1,size(Gs,3))), 1, [])';

%% Case 5: Small Signal Dependent Noise, Long Head Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0   ,0; 0,0.01];

% The length of the holding period
Thold = 0.100 / Delta;
figNum = 5;
figName = 'Small Signal Dependent Noise, Long Head Hold';
% Perform the simulations and plot the results
[Ks, Gs, eye_mean, head_mean, gaze_mean] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, figNum, figName);

kalmanGains(:,5) = reshape(mat2cell(Ks, size(Ks,1), size(Ks,2), ones(1,size(Ks,3))), 1, [])';
feedbackGains(:,5) = reshape(mat2cell(Gs, size(Gs,1), size(Gs,2), ones(1,size(Gs,3))), 1, [])';
