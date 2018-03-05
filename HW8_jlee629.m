function HW8_jlee629()


%% part 1 Fourier Domain filtering
X = imread('Nice','jpg');
X = rgb2gray(X);
Y = fft2(X);
%FFT of the initial image
Z1 = fftshift(Y);
Z1 = abs(Z1);
Z1 = log(Z1+1);
Z1 = mat2gray(Z1);

[row, col] = size(Y);

%High pass filter. 
highpass = 0.1; % value from 0 to 1, with 0 implying no filter, 0.1 indicates that the lower 10% of the frequencies are blocked.
Z2 = Y;

Z2 = fftshift(Z2);
Z2(floor(row*(0.5-highpass/2)):floor(row*(0.5+highpass/2)),floor(col*(0.5-highpass/2)):floor(col*(0.5+highpass/2)))=0;
X2 = ifftshift(Z2);
X2 = ifft2(X2);
X2 = abs(X2);
Z2 = log(abs(Z2)+1);
Z2 = mat2gray(Z2);
Z2 = abs(Z2);



%Low pass filter
Lowpass = 0.9;
Z3 = Y;

Z3(floor(row*(0.5-Lowpass/2)):floor(row*(0.5+Lowpass/2)),:)=0;
Z3(:,floor(col*(0.5-Lowpass/2)):floor(col*(0.5+Lowpass/2)))=0;
X3 = ifft2(Z3);
X3 = abs(X3);
Z3 = fftshift(Z3);
Z3 = log(abs(Z3)+1);
Z3 = mat2gray(Z3);
Z3 = abs(Z3);


% Bandpass filter

Bandpass = [0.4 0.5];

Z4 = Y;
Z4 = fftshift(Z4);
for n = 1:row
    for m = 1:col
        if m >= -Bandpass(2)*col/2 + col/2 && m <= Bandpass(2)*col/2 + col/2
            if n>abs(Bandpass(1)*row/2+row/2) && n<abs(Bandpass(2)*row/2+row/2)               
                Z4(n,m) = 0;
            elseif n<-Bandpass(1)*row/2+row/2 && n>-Bandpass(2)*row/2+row/2
                Z4(n,m) = 0;
            end
        end
    end
end



for m = 1:col
    for n = 1:row
        if n >= -Bandpass(1)*row/2 + row/2 && n <= Bandpass(2)*row/2 + row/2
            if m>abs(Bandpass(1)*col/2+col/2) && m<abs(Bandpass(2)*col/2+col/2)
                Z4(n,m) = 0;
            elseif m<-Bandpass(1)*col/2+col/2 && m>-Bandpass(2)*col/2+col/2
                Z4(n,m) = 0;
            end
        end
        
    end
end


X4 = abs(ifft2(ifftshift(Z4)));
X4 = mat2gray(X4);

Z4 = log(abs(Z4)+1);
Z4 = abs(Z4);
Z4 = mat2gray(Z4);




%plotting figures

figure
colormap('gray');
subplot(2,2,1)
imagesc(X)
title('real image')
subplot(2,2,2)
imagesc(X2)
title('High pass filter')
subplot(2,2,3)
imagesc(X3)
title('Low pass filter')
subplot(2,2,4)
imagesc(X4)
title('Band pass filter')


figure
colormap('gray');
subplot(2,2,1)
imagesc(Z1)
title('real image')
subplot(2,2,2)
imagesc(Z2)
title('High pass filter')
subplot(2,2,3)
imagesc(Z3)
title('Low pass filter')
subplot(2,2,4)
imagesc(Z4)
title('Band pass filter')


%% part 2 Real and Imaginary parts of the DFT. 


X1 = imread('foret','jpg');
X2 = imread('Nice','jpg');
X1 = rgb2gray(X1);
X2 = rgb2gray(X2);


Y1 = fft2(X1);
Y2 = fft2(X2);

%extracting magnitude and phase of each fourier transform
Z.Y1mag = abs(fftshift(Y1));
Z.Y1phase = angle(fftshift(Y1));
Z.Y2mag = abs(fftshift(Y2));
Z.Y2phase = angle(fftshift(Y2));


%swapping magnitude and phase of each fourier transform
Mod1 = Z.Y1mag.*exp(1i*Z.Y2phase) ;
Mod2 = Z.Y2mag.*exp(1i*Z.Y1phase) ;
%inverse fourier transformation
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


%% homework 7-1

dt = 0.1*1e-3;
alpha = 0.01;
t = [dt : dt :1];
for n = 1:length(t)
    h(n) = exp(-alpha*t(n));
end
    


plot(h)

















