function HW3_jlee629()


Dataset = load('classify_regress.dat');

N = length(Dataset);
%initializing parameters

X = [ones(N,1) Dataset(:,1:2)];

%Calculating the weight vector using the normal equation.
w = (X.'*X)^-1*X.'*Dataset(:,3);
disp(w)

%plot Dataset. red if in class 0, green if in class 1
figure

scatter(Dataset(find(Dataset(:,3)==0),1),Dataset(find(Dataset(:,3)==0),2),'r')
hold on
scatter(Dataset(find(Dataset(:,3)==1),1),Dataset(find(Dataset(:,3)==1),2),'g')

xlabel('x1')
ylabel('x2')

xaxis = [-6:1:10];
for x = 1:length(xaxis)
yaxis(x) = -1*(w(1)+w(2)*xaxis(x)-0.5)/w(3);
end

plot(xaxis,yaxis);
legend('class 0', 'class 1','classification boundary');