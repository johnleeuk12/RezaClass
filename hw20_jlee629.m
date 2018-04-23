function hw20_jlee629()

%% Initializing Parameters 
clear all
Data = load('brainimage.txt');
Data = Data/max(max(Data));

figure
imshow(Data)
title('Image of a brain')

Data = reshape(Data,[],1);
Sigma = std(Data);
mu = [0.3 0.5 0.7];
P_z = [1./3 1./3 1./3]; % initial probabilities for each class
Nb_class = 3;


P_x_knowing_z = zeros(length(Data),Nb_class);

% for step = 1:100
% E step
    for n = 1:length(Data)
        for i = 1:Nb_class
            P_x_knowing_z(n,i) = sqrt(2*pi*Sigma)^-1*exp(-0.5*(Data(n)-mu(i))^2/Sigma);            
        end
        for i = 1:Nb_class
        P_z_knowing_x(n,i) = P_x_knowing_z(n,i)*P_z(i)/(sum(P_x_knowing_z(n,:).*P_z)); % calculate the posterior probability
        end
    end
        
% end