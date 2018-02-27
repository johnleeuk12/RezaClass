function hw5_jlee629()


%generate dataset 
X = rand(2500,1)*10;
X = [ones(2500,1) X X.*X];

w = [0.2 0.5 0.1];

y = X*w.';
noise = randn(2500,1);

% we observe y(i) + noise(i) for i = 1:2500



meansigmaML = [];
for N = 5:1000
    if mod(N,100)==0
    disp(N)
    end
    sigmaML = [];
    for i = 1:100 % trial repetitions
        indx = randsample([1:2500],N);
        
        sigmaML = [sigmaML mean(noise(indx).^2)];
    end
    
    meansigmaML  = [meansigmaML mean(sigmaML)];
end


figure

plot(meansigmaML)