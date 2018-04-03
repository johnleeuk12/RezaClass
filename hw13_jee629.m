function hw13_jee629()


x = 0:0.1:300;
px_bass = 0.6*pdf('Normal',x,165,7)+0.4*pdf('Normal',x,180,7);
px_salm = 0.6*pdf('Normal',x,180,3)+0.4*pdf('Normal',x,160,8);

figure
plot(x,px_bass)
hold on
plot(x,px_salm)
axis([130, 210, 0, inf]);
xlabel('luminosity')
title('Likelihood')
legend('Bass','Salmon')



% marginal distribution
px = px_bass*0.4+px_salm*0.6;

figure
plot(x,px)
axis([130, 210, 0, inf]);
xlabel('luminosity')
title('Marginal distribution')
% legend('Bass','Salmon')


% posterior probabilities for each fish species
pc_bass = px_bass*0.4./px;
pc_salm = px_salm*0.6./px;

figure

plot(x,pc_bass)
hold on
plot(x,pc_salm)
axis([130, 210, 0, inf]);
xlabel('luminosity')
title('posterior probabilities')
legend('Bass','Salmon')


%Calculating the variance of the posterior probability distribution

var = pc_bass.*pc_salm;

figure
plot(x,var)
axis([130, 210, 0, inf]);
xlabel('luminosity')
title('Variance of distribution')




%% Homework 14
clear all
Data = load('classification.dat');


Variance = zeros(2,2);
for cls = 1:3
    N(cls) = length(find(Data(:,1) == cls));
    
    mu(cls,:) = mean(Data(find(Data(:,1) == cls),2:3));
    
    Sigma{cls} = zeros(2,2);
    for n = find(Data(:,1) == cls)
        Sigma{cls} = Sigma{cls}+(Data(n,2:3)-mu(cls,:)).'*(Data(n,2:3)-mu(cls,:));
    end
    Variance = Variance+Sigma{cls};
    Sigma{cls} = Sigma{cls}/N(cls);
end    
Variance = Variance/length(Data);


% for boundary 1/2
W(1,:) = mu(1,:)*inv(Variance)-mu(2,:)*inv(Variance);
W(2,:) = mu(2,:)*inv(Variance)-mu(3,:)*inv(Variance);
W(3,:) = mu(3,:)*inv(Variance)-mu(1,:)*inv(Variance);
Wo(1) = 0.5*(-mu(1,:)*inv(Variance)*mu(1,:).' + mu(2,:)*inv(Variance)*mu(2,:).');
Wo(2) = 0.5*(-mu(2,:)*inv(Variance)*mu(2,:).' + mu(3,:)*inv(Variance)*mu(3,:).');   
Wo(3) = 0.5*(-mu(3,:)*inv(Variance)*mu(3,:).' + mu(1,:)*inv(Variance)*mu(1,:).');   
    
x2 = -10:0.5:10;

figure

for cls = 1:3
    x1{cls} = -Wo(cls)/W(cls,2)-W(cls,1)/W(cls,2)*x2;
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
    plot(x2,x1{cls});
    axis([-8,2,-8,6]);
end
vec = -10:0.1:10;
PClass = zeros(length(vec),length(vec));

for n1 = 1: length(vec)
    for n2 = 1: length(vec)
        xx1 = vec(n1);
        xx2 = vec(n2);
        P(1) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(1,:))*inv(Variance)*([xx1,xx2]-mu(1,:)).'));
        P(2) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(2,:))*inv(Variance)*([xx1,xx2]-mu(2,:)).'));    
        P(3) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(3,:))*inv(Variance)*([xx1,xx2]-mu(3,:)).'));
        PClass(n1,n2) = find(P==max(P));
    end
end


%Quadratic

WWq{1}= -0.5*(inv(Sigma{1})-inv(Sigma{2}));   
WWq{2}= -0.5*(inv(Sigma{2})-inv(Sigma{3}));      
WWq{3}= -0.5*(inv(Sigma{3})-inv(Sigma{1}));       

Wq(1,:) = mu(1,:)*inv(Sigma{1})-mu(2,:)*inv(Sigma{2});
Wq(2,:) = mu(2,:)*inv(Sigma{2})-mu(3,:)*inv(Sigma{3});
Wq(3,:) = mu(3,:)*inv(Sigma{3})-mu(1,:)*inv(Sigma{1});

Woq(1) = 0.5*(-mu(1,:)*inv(Sigma{1})*mu(1,:).' + mu(2,:)*inv(Sigma{2})*mu(2,:).') -0.5*log(norm(Sigma{1})/norm(Sigma{2}));
Woq(2) = 0.5*(-mu(2,:)*inv(Sigma{2})*mu(2,:).' + mu(3,:)*inv(Sigma{3})*mu(3,:).') -0.5*log(norm(Sigma{2})/norm(Sigma{3}));
Woq(3) = 0.5*(-mu(3,:)*inv(Sigma{3})*mu(3,:).' + mu(1,:)*inv(Sigma{1})*mu(1,:).') -0.5*log(norm(Sigma{3})/norm(Sigma{1}));    


syms x1 x2

eqn{1} = [x1 x2]*WWq{1}*[x1 x2].' + Wq(1,:)*[x1 x2].' + Woq(1);
eqn{2} = [x1 x2]*WWq{2}*[x1 x2].' + Wq(2,:)*[x1 x2].' + Woq(2);
eqn{3} = [x1 x2]*WWq{3}*[x1 x2].' + Wq(3,:)*[x1 x2].' + Woq(3);
    

figure

for cls = 1:3
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
    fimplicit(eqn{cls})
    axis([-8,2,-8,6]);
end   
    
    
% Gaussian Kernel density estimation 

syms u
KK(u) = 1/2*pi*exp(-0.5*u.'*u);
    
eqnGaus{1} = 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





