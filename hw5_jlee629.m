function hw5_jlee629()


%Initializing parameters and dataset 

X = rand(2500,1)*10;
X = [ones(2500,1) X X.*X].';
W = [0.2 0.5 0.1].';

Y = W.'*X;
noise = randn(1,2500);
sigmaML = zeros(1,length(5:500));
nn = 1 ;



%  running simulation to calculate the sigma to obtain the maximum
%  likelyhood
for N = 5:500 % N is the number of sampling points
    if mod(N,50)==0
        disp(N)
    end
    for i = 1:100 % for each N, a 100 iterations
    indx = randsample([1:2500],N);
    sigmaML(nn) = sigmaML(nn)+mean(noise(indx).^2);
    end
    sigmaML(nn) = sigmaML(nn)/100;
    nn = nn+1;
end




figure
title('\sigma^2_{ML} depending on the number of sampling points') 
ylabel('\sigma^2_{ML}')
xlabel('The number of sampling points N')
plot(sigmaML)