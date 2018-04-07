function hw16_jlee629()

% Inputs

X = [[0 1 0 1]' [0 0 1 1]'];
y = [0 1 1 0]';

param = rand(1,7);
sigmo = @(x) 1./(1+exp(-x));
eta = 0.2;

%estimating parameters

var1 = param(4)*X(:,1)+param(6) + param(5)*X(:,2)+param(6);
var2 = param(1)*X(:,1)+param(6) + param(2)*X(:,2)+param(6) + param(3)*sigmo(var1);
yhat = sigmo(var2);
Err = y-yhat;
Dparam = zeros(1,7);
for p = 1:4
    Dparam(1) = Dparam(1)+ Err(p)*yhat(1-yhat)*X(p,1);
    Dparam(2) = Dparam(2)+ Err(p)*yhat(1-yhat)*X(p,2);
    Dparam(3) = Dparam(3)+ Err(p)*yhat(1-yhat)*sigmo(var1(p));
    Dparam(4) = Dparam(4)+ Err(p)*yhat(1-yhat)*
    
    
    http://bengio.abracadoudou.com/lectures/mlp.pdf     
    http://bengio.abracadoudou.com/lectures/old/tex_ann.pdf
    

