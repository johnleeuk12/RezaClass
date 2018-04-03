function hw14_jlee629()


clear all
data = load('classification.dat');


Variance = zeros(2,2);
for cls = 1:3
    % Calculate the number of data for each class
    N(1,cls) = length(find(data(:,1) == cls));
    
    % Calculate the mean
    mu(cls,:) = [sum(data(find(data(:,1) == cls),2)) sum(data(find(data(:,1) == cls),3))];
    mu(cls,:) = mu(cls,:)/N(1,cls);
    
    % Calculate the standard deviation.
    Sigma{cls} = zeros(2,2);
    for n = find(data(:,1) == cls)
        Sigma{cls} = Sigma{cls} + (data(n,2:3)-mu(cls,:)).'*(data(n,2:3)-mu(cls,:));
    end
    Sigma{cls} = Sigma{cls}/N(1,cls);    
    
    %Calculate the single variance
    Variance = Variance + N(1,cls)*Sigma{cls};
end

Variance = Variance/sum(N);

X = data(:,2:3);

W = Variance*mu(1,:).' + Variance*mu(2,:).' ;
W0 = -0.5*mu(1,:)*Variance*mu(1,:).'+-0.5*mu(2,:)*Variance*mu(2,:).';

Y = X*W+W0;