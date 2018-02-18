function HW3()

clear all

Dataset = load('classify_regress.dat'); %{x1,x2,y}
N = size(Dataset,1);
%using steepest descent learning rule for data classification.

%% initializing variables
w = [rand rand rand];
accuracy = 0;
eta = 0.5;


while accuracy(end) < 0.9
    delta = [0 0 0];
   error = 0;
    for n = 1:N
         if w*[1 Dataset(n,1:2)].'> 0.5
             yhat = 1;
         else
             yhat = 0;
         end
         if yhat ~= Dataset(n,3)
             error = error+1;
         end
        delta = delta + [1 Dataset(n,1:2)]*(Dataset(n,3)-yhat);
    end
    w = w + eta*1/size(Dataset,1)*delta;
    accuracy = [accuracy (N-error)/N];
end
    


figure
scatter(Dataset([find(Dataset(:,3)==0)],1),Dataset([find(Dataset(:,3)==0)],2),[],'r')
hold on
scatter(Dataset([find(Dataset(:,3)==1)],1),Dataset([find(Dataset(:,3)==1)],2),[],'g')