function hw14_jlee629()

clear all

Data = load('classification.dat');

%Calculate the mean and variance of each class
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


% Take a number of points in the data space, and classify each point.
vec = -10:0.1:10;
PClass = zeros(length(vec),length(vec));

for n1 = 1: length(vec)
    for n2 = 1: length(vec)
        xx1 = vec(n1);
        xx2 = vec(n2);
        P(1) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(1,:))*inv(Variance)*([xx1,xx2]-mu(1,:)).'));
        P(2) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(2,:))*inv(Variance)*([xx1,xx2]-mu(2,:)).'));    
        P(3) = 1/sqrt(norm(Variance))*exp(-.5*(([xx1,xx2]-mu(3,:))*inv(Variance)*([xx1,xx2]-mu(3,:)).'));
        PClass(n1,n2) = find(P==max(P))*0.2+0.5;
    end
end
figure
imagesc(vec,vec,PClass.')
axis([-8,2,-8,6]);
colormap(gray)
caxis([0 1])
set(gca,'YDir','normal')
hold on


% Calculate the boundary between each class

W(1,:) = mu(1,:)*inv(Variance)-mu(2,:)*inv(Variance);
W(2,:) = mu(2,:)*inv(Variance)-mu(3,:)*inv(Variance);
W(3,:) = mu(3,:)*inv(Variance)-mu(1,:)*inv(Variance);
Wo(1) = 0.5*(-mu(1,:)*inv(Variance)*mu(1,:).' + mu(2,:)*inv(Variance)*mu(2,:).');
Wo(2) = 0.5*(-mu(2,:)*inv(Variance)*mu(2,:).' + mu(3,:)*inv(Variance)*mu(3,:).');   
Wo(3) = 0.5*(-mu(3,:)*inv(Variance)*mu(3,:).' + mu(1,:)*inv(Variance)*mu(1,:).');  

x2 = -10:0.5:10;
for cls = 1:3
    x1{cls} = -Wo(cls)/W(cls,2)-W(cls,1)/W(cls,2)*x2;
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
    plot(x2,x1{cls});
    axis([-8,2,-8,6]);
end

hold off


%% Quadratic

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

PclassQuad = zeros(length(vec),length(vec));
for n1 = 1: length(vec)
    for n2 = 1: length(vec)
        xx1 = vec(n1);
        xx2 = vec(n2);
        P(1) = 1/sqrt(norm(Sigma{1}))*exp(-.5*(([xx1,xx2]-mu(1,:))*inv(Sigma{1})*([xx1,xx2]-mu(1,:)).'));
        P(2) = 1/sqrt(norm(Sigma{2}))*exp(-.5*(([xx1,xx2]-mu(2,:))*inv(Sigma{2})*([xx1,xx2]-mu(2,:)).'));    
        P(3) = 1/sqrt(norm(Sigma{3}))*exp(-.5*(([xx1,xx2]-mu(3,:))*inv(Sigma{3})*([xx1,xx2]-mu(3,:)).'));
        PClassQuad(n1,n2) = find(P==max(P))*0.2+0.5;
    end
end


figure
imagesc(vec,vec,PClassQuad.')
axis([-8,2,-8,6]);
colormap(gray)
caxis([0 1])
hold on
set(gca,'YDir','normal')

for cls = 1:3
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
    fimplicit(eqn{cls})
    axis([-8,2,-8,6]);
end   
    
    
%% Gaussian Kernel density estimation 

h = .5;
PclassGauss = zeros(length(vec),length(vec));

for n1 = 1: length(vec)
    for n2 = 1: length(vec)
        xx1 = vec(n1);
        xx2 = vec(n2);
        P = zeros(1,3);
        for cls = 1:3
            ind = find(Data(:,1) == cls);
            for nn = 1:length(ind)
                u = ([xx1 xx2]-Data(ind(nn),2:3))/h;
                P(cls) = P(cls)+1/(2*pi)*exp(-0.5*(u*u.'));
            end
        end
%         P(2) = 1/sqrt(norm(Sigma{2}))*exp(-.5*(([xx1,xx2]-mu(2,:))*inv(Sigma{2})*([xx1,xx2]-mu(2,:)).'));    
%         P(3) = 1/sqrt(norm(Sigma{3}))*exp(-.5*(([xx1,xx2]-mu(3,:))*inv(Sigma{3})*([xx1,xx2]-mu(3,:)).'));
        PclassGauss(n1,n2) = find(P==max(P))*0.2+0.5;
    end
end

 
figure
imagesc(vec,vec,PclassGauss.')
axis([-8,2,-8,6]);
colormap(gray)
caxis([0 1])
hold on
set(gca,'YDir','normal')

for cls = 1:3
    scatter(Data(find(Data(:,1) == cls),2),Data(find(Data(:,1) == cls),3),20,'filled')
    hold on
    axis([-8,2,-8,6]);
end   
        
    
    
    
    
    
    
    
    
    
    
    
    
    
