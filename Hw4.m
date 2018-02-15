function Hw4()

%initializing parameters

clear all

c = 0;
global sigma
sigma = 0.7;
x = [-10:0.1:10];
w = 3;
g =  exp(-1/(2*sigma^2)*(x-c*ones(1,length(x))).^2);

figure
title('Gaussian bases')
plot(x,w*g)
hold on 
plot(x,g)

xlabel('x')
ylabel('g')
legend('weight = 3', 'weight = 1')
Dataset = [1 1 ; -1 1];


J = [];
for c = x
    x1 = Dataset(1,1);
    x2 = Dataset(2,1);
    y1hat =  w*gaussi(x1,c);
    y2hat =  w*gaussi(x2,c);
    J = [J 1/2*((Dataset(1,2)-y1hat)^2+(Dataset(2,2)-y2hat)^2)];
end





figure 
plot(x,J)
hold on 
title('Loss function for c')
ylabel('squared error lost function J')
xlabel('x')  


w1 = 1;
w2 = 1;

J2 =[];
for c = x
    x1 = Dataset(1,1);
    x2 = Dataset(2,1);
    y1hat1 =  w1*gaussi(x1,1);
    y2hat1 =  w1*gaussi(x2,1);
    y1hat2 =  w2*gaussi(x1,c);
    y2hat2 =  w2*gaussi(x2,c);   
    J2 = [J2 1/2*((Dataset(1,2)-y1hat1-y1hat2)^2 ...
    +(Dataset(2,2)-y2hat1-y2hat2)^2)];
end

figure 
plot(x,J2)
hold on 
title('Loss function for c2')
ylabel('squared error lost function J')
xlabel('x')  





function out = gaussi(x,c)
global sigma
out = exp(-1/(2*sigma^2)*(x-c*ones(1,length(x))).^2);
