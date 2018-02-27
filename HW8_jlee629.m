function HW8_jlee629()

X = imread('Paris','jpg');
X = rgb2gray(X);
figure
imagesc(X)
colormap('gray');

Y = fft2(X);
Z1 = fftshift(Y);
X1 = ifft2(Y);
figure
colormap('gray');
imagesc(X1);

Z1 = abs(Z1);
Z1 = log(Z1+1);
Z1 = mat2gray(Z1);
figure
imshow(Z1)

Lowpass = 0.1; % value from 0 to 1, with 0 implying no filter. 
Z2 = Y;
% Z2(1:Lowpass*size(Y,1),1:Lowpass*size(Y,2)) = 0;
Z2 = fftshift(Z2);
Z2(floor(size(Y,1)*(0.5-Lowpass/2)):floor(size(Y,1)*(0.5+Lowpass/2)),floor(size(Y,2)*(0.5-Lowpass/2)):floor(size(Y,2)*(0.5+Lowpass/2)))=0;
X2 = ifftshift(Z2);
X2 = ifft2(X2);
X2 = abs(X2);
Z2 = log(abs(Z2)+1);
Z2 = mat2gray(Z2);
Z2 = abs(Z2);
figure
imshow(Z2)

figure
colormap('gray');
imagesc(X2)

Highpass = 0.9;
Z3 = Y;
% Z2(1:Lowpass*size(Y,1),1:Lowpass*size(Y,2)) = 0;


Z3(floor(size(Y,1)*(0.5-Highpass/2)):floor(size(Y,1)*(0.5+Highpass/2)),:)=0;
Z3(:,floor(size(Y,2)*(0.5-Highpass/2)):floor(size(Y,2)*(0.5+Highpass/2)))=0;
X3 = ifft2(Z3);
X3 = abs(X3);
Z3 = fftshift(Z3);
Z3 = log(abs(Z3)+1);
Z3 = mat2gray(Z3);
Z3 = abs(Z3);
figure
imshow(Z3)

figure 
colormap('gray')
subplot(1,3,1)
imagesc(X)
subplot(1,3,2)
imagesc(X2)
subplot(1,3,3)
imagesc(X3)