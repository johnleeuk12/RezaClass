function hw19_jlee629()

%% Parameters
N = 64;
X = [];
m = 1; % individual number
for n = 1:N
    img = loadimage(m,n);
    X = [X reshape(img,[],1)];
end
D = length(X);

%% apply PCA
d = 2;
[mu,Ud,Y] = pcaJL(X,d);
img_average = reshape(mu,[192,168]);
eface1 = reshape(Ud(:,1),192,168);
eface2 = reshape(Ud(:,2),192,168);

figure 
imshow(uint8(img_average))
figure
imshow(eface1,[min(eface1(:)),max(eface1(:))])
figure
imshow(eface2,[min(eface2(:)),max(eface2(:))])


y1=-std(Y(1,:)):0.2*std(Y(1,:)):std(Y(1,:));
Z1=mu+Ud(:,1)*y1;
y2=-std(Y(2,:)):0.2*std(Y(2,:)):std(Y(2,:));
Z2=mu+Ud(:,2)*y2;

figure
for i = 1:length(y1)
    subplot(3,4,i)
    imshow(uint8(reshape(Z1(:,i),192,168)))
end

figure
for i = 1:length(y2)
    subplot(3,4,i)
    imshow(uint8(reshape(Z2(:,i),192,168)))
end




%% PCA function
function [mu,Ud,Y] = pcaJL(X,d)

mu = mean(X,2);
centerX = X-mu;
[U,~,~] = svd(centerX,'econ');
Ud=U(:,1:d);
Y=Ud'*centerX;



% img2 = reshape(mu,[192,168]);