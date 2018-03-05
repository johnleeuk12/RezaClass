function HW8_jlee629()


%% part 1 Fourier Domain filtering
X = imread('Nice','jpg');
X = rgb2gray(X);
figure
imagesc(X)
colormap('gray');

% Y = fft2(X);
% Z1 = fftshift(Y);
% X1 = ifft2(Y);
% figure
% colormap('gray');
% imagesc(X1);

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


Bandpass = [0.2 0.4];
Z4 = Y;


Z4 = fftshift(Z4);
for x = 1:size(Z4,1)
    for y = 1:size(Z4,2)
        if ((x-400)^2 + (y-640)^2)^0.5 > Bandpass(1)*size(Z4,1) && ((x-400)^2 + (y-640)^2)^0.5 < Bandpass(2)*size(Z4,1)
            Z4(x,y) = 0;
        end
    end
end
X4 = ifftshift(Z4);
X4 = ifft2(X4);
X4 = abs(X4);
Z4 = log(abs(Z4)+1);
Z4 = mat2gray(Z4);
Z4= abs(Z4);
figure
imshow(Z4)




figure 
colormap('gray')
subplot(2,2,1)
imagesc(X)
subplot(2,2,2)
imagesc(X2)
subplot(2,2,3)
imagesc(X3)
subplot(2,2,4)
imagesc(X4)


%% part 2 Real and Imaginary parts of the DFT. 


X1 = imread('foret','jpg');
X2 = imread('Nice','jpg');
X1 = rgb2gray(X1);
X2 = rgb2gray(X2);


Y1 = fft2(X1);
Y2 = fft2(X2);

Z.Y1mag = abs(fftshift(Y1));
Z.Y1phase = angle(fftshift(Y1));
Z.Y2mag = abs(fftshift(Y2));
Z.Y2phase = angle(fftshift(Y2));



Mod1 = Z.Y1mag.*exp(1i*Z.Y2phase) ;
Mod2 = Z.Y2mag.*exp(1i*Z.Y1phase) ;
XMod1 = abs(ifft2(Mod1));
XMod2 = abs(ifft2(Mod2));
fourier1 = log(abs(fftshift(Mod1))+1);
fourier2 = log(abs(fftshift(Mod2))+1);

figure
subplot(2,2,1)
imshow(abs(mat2gray(log(abs(Z.Y1mag+1)))))
title('Magnitude(Forest)')
subplot(2,2,2)
imshow(abs(mat2gray(log(abs(Z.Y2mag+1)))))
title('Magnitude(Room)')
subplot(2,2,3)
imshow(abs(mat2gray(log(abs(Z.Y1phase+1)))))
title('Phase(Forest)')
subplot(2,2,4)
imshow(abs(mat2gray(log(abs(Z.Y2phase+1)))))
title('Phase(Room)')


figure
colormap('gray')
subplot(2,2,1)

imagesc(X1)
title('Forest')
subplot(2,2,2)

imagesc(X2)
title('Room')
subplot(2,2,3)

imagesc(XMod1)
title('Magnitude(Forest)+Phase(Room)')
subplot(2,2,4)
imagesc(XMod2)
title('Magnitude(Room)+Phase(Forest)')






















