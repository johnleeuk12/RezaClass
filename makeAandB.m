function [A,B] = makeAandB(tau1eye, tau2eye, alpha1eye, tau1head,tau2head, alpha1head , Delta)
%%
    Ae = [0 1 0; ...
        -(tau1eye*tau2eye)^-1 -(tau1eye+tau2eye)*(tau1eye*tau2eye)^-1 (tau1eye*tau2eye)^-1; ...
        0 0 -(alpha1eye)^-1];
    be = [0 0 alpha1eye^-1].';
    Ah = [0 1 0; ...
        -(tau1head*tau2head)^-1 -(tau1head+tau2head)*(tau1head*tau2head)^-1 (tau1head*tau2head)^-1; ...
        0 0 -(alpha1head)^-1];
    bh = [0 0 alpha1head^-1].';

    Ac = [Ae zeros(3,3);zeros(3,3) Ah];
    Bc = [be zeros(3,1); zeros(3,1) bh];


    A = [expm(Ac*Delta) zeros(6,1); zeros(1,6) 1];
    B = [Ac\(expm(Ac*Delta)-eye(size(Ac)))*Bc ; zeros(1,2)];

end