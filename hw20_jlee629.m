function hw20_jlee629()

%% Initializing Parameters 
clear all
Data = load('brainimage.txt');
Data = Data/max(max(Data));
% 
% figure
% imshow(Data)
% title('Image of a brain')
[x1, x2] = size(Data);
Data = reshape(Data,[],1);
Sigma_hat = std(Data)*ones(1,3);
mu_hat = [0.3 0.5 0.7];
P_z = [1./3 1./3 1./3]; % initial probabilities for each class
Nb_class = 3;
Likelyhood = [];

P_x_knowing_z = zeros(length(Data),Nb_class);
% tau = zeros(length(Data),Nb_class);
n_hat = zeros(1,Nb_class);

%% EM estimation 
for step = 1:50
    % E step, estimating probabilities. 
    for n = 1:length(Data)
        for i = 1:Nb_class
            P_x_knowing_z(n,i) = sqrt(2*pi*Sigma_hat(i))^-1*exp(-0.5*(Data(n)-mu_hat(i))^2/Sigma_hat(i));
        end
        for i = 1:Nb_class
            P_z_knowing_x(n,i) = P_x_knowing_z(n,i)*P_z(i)/(sum(P_x_knowing_z(n,:).*P_z)); % calculate the posterior probability
        end
    end
    Likelyhood = [Likelyhood sum(log(P_x_knowing_z*P_z.'))];
    % M step, Re-estimating parameters. 
    for i = 1:Nb_class
        n_hat(i) = sum(P_z_knowing_x(:,i));
        P_z(i) = n_hat(i)/length(Data);
        mu_hat(i) = 1./n_hat(i)*sum(P_z_knowing_x(:,i).*Data);
        Sigma_hat(i) = 1./n_hat(i)*sum(P_z_knowing_x(:,i).*(Data-mu_hat(i)).^2);
    end
end


figure
plot(Likelyhood)
xlabel('steps')
title('Likelyhood')

% Classifying each pixel based on the posterior probability.
Class_estimation = zeros(length(Data),1);
Class1 = zeros(length(Data),1);
Class2 = zeros(length(Data),1);
Class3 = zeros(length(Data),1);

for n = 1:length(Data)
    Class_estimation(n,1) = find(P_z_knowing_x(n,:) == max(P_z_knowing_x(n,:)));
    %plotting each class
    if Class_estimation(n,1) == 1
        Class1(n,1) = 1;
    elseif Class_estimation(n,1) == 2
        Class2(n,1) = 2;
    else
        Class3(n,1) = 3;
    end
end

Class_estimation = reshape(Class_estimation,x1,x2);
Class1 = reshape(Class1,x1,x2);
Class2 = reshape(Class2,x1,x2);
Class3 = reshape(Class3,x1,x2);

Class_estimation = Class_estimation/max(max(Class_estimation)); 
Class1 = Class1/max(max(Class1));
Class2 = Class2/max(max(Class2));
Class3 = Class3/max(max(Class3));

figure

subplot(2,2,1)
imshow(Class_estimation)
xlabel('reconstructed image')

subplot(2,2,2)
imshow(Class1)
xlabel('Class 1')

subplot(2,2,3)
imshow(Class2)
xlabel('Class 2')

subplot(2,2,4)
imshow(Class3)
xlabel('Class 3')

