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









