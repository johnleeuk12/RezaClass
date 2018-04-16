function hw19_jlee629()

%% Parameters
N = 64;
X = [];
for n = 1:N
    img = loadimage(1,n);
    X = [X reshape(img,[],1)];
end
D = length(X);

%% apply PCA
d = 2;
[mu,Ud,Y] = pcaJL(X,d);
img_average = reshape(mu,[192,168]);
eface1 = 

%% PCA function
function [mu,Ud,Y] = pcaJL(X,d)

mu = mean(X,2);
centerX = X-mu;
[U,~,~] = svd(centerX);
Ud=U(:,1:d);
Y=Ud'*CenterX;



% img2 = reshape(mu,[192,168]);