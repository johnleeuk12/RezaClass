function HW2SB2()

Dataset = load('ImageData.mat');
X = Dataset.X;
Y = Dataset.Y;



%% Optimal linear mapping

%Calculating W and the rank of W.
W = Y*pinv(X);
k = rank(W);


% estimating Y using W and input matrix X 
Yhat = W*X;

figure
colormap(gray)
for i = 1:10
    subplot(2,5,i)
    m = reshape(Yhat(:,i),33,33);
    imagesc(m)
    hold on
end

hold off

% corrupting input matrix X

Xdash = X;
l = 112; %length of each input image.
h = 80;  %height of each input image.

for i = 1:10
    for n = 1:h
        Xdash(56*h+1:end,i) = 0;
    end
end
figure
colormap(gray)

for i = 1:10
    subplot(2,5,i)
m=reshape(Xdash(:,i), 80, 112);
imagesc(m)
end

Yhatnew = W*Xdash;

figure
colormap(gray)
for i = 1:10
    subplot(2,5,i)
    m = reshape(Yhatnew(:,i),33,33);
    imagesc(m)
    hold on
end

hold off

% Corrupting weight matrix W
W(randsample(1:length(X)*length(Y),length(X)*length(Y)*0.9)) = 0;
Yhat = W*X;

figure
colormap(gray)
for i = 1:10
    subplot(2,5,i)
    m = reshape(Yhat(:,i),33,33);
    imagesc(m)
    hold on
end

hold off




%% Auto associative memory

Wa = X*pinv(X);
ka = rank(W);

%Corrupting X

Yhata = Wa*Xdash;

figure
colormap(gray)
for i = 1:10
    subplot(2,5,i)
    m = reshape(Yhata(:,i),80, 112);
    imagesc(m)
    hold on
end

% Corrupting weight matrix W
Wa(randsample(1:length(X)*length(Y),length(X)*length(Y)*0.9)) = 0;
Yhata = Wa*X;

figure
colormap(gray)
for i = 1:10
    subplot(2,5,i)
    m = reshape(Yhata(:,i),80, 112);
    imagesc(m)
    hold on
end

hold off





