function hw16_jlee629()

clear all
% Inputs

X = [[0 1 0 1]' [0 0 1 1]'];
y = [0 1 1 0]';

param = rand(1,7)*10;
sigmo = @(x) 1./(1+exp(-x));
eta = 2;

%estimating parameters
for trial = 1:100000
    var1 = param(4)*X(:,1) + param(5)*X(:,2)+param(7);
    var2 = param(1)*X(:,1) + param(2)*X(:,2)+param(6) + param(3)*sigmo(var1);
    yhat = sigmo(var2);
    Err = y-yhat;
    LMS(trial) = immse(y,yhat);
    Dparam = zeros(1,7);
    for p = 1:4
        Dparam(1) = Dparam(1)+ Err(p)*yhat(p)*(1-yhat(p))*X(p,1);
        Dparam(2) = Dparam(2)+ Err(p)*yhat(p)*(1-yhat(p))*X(p,2);
        Dparam(3) = Dparam(3)+ Err(p)*yhat(p)*(1-yhat(p))*sigmo(var1(p));
        Dparam(4) = Dparam(4)+ Err(p)*yhat(p)*(1-yhat(p))*param(3)*sigmo(var1(p))*(1-sigmo(var1(p)))*X(p,1);
        Dparam(5) = Dparam(5)+ Err(p)*yhat(p)*(1-yhat(p))*param(3)*sigmo(var1(p))*(1-sigmo(var1(p)))*X(p,2);
        Dparam(6) = Dparam(6)+ Err(p)*yhat(p)*(1-yhat(p));
        Dparam(7) = Dparam(7)+ Err(p)*yhat(p)*(1-yhat(p))*param(3)*sigmo(var1(p))*(1-sigmo(var1(p)));
    end
    param = param + eta*Dparam;
end
    

plot(LMS)