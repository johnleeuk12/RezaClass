function hw13_jee629()


x = 0:0.1:300;
px_bass = 0.6*pdf('Normal',x,165,7)+0.4*pdf('Normal',x,180,7);
px_salm = 0.6*pdf('Normal',x,180,3)+0.4*pdf('Normal',x,160,8);

plot(x,px_bass)
hold on
plot(x,px_salm)

% marginal distribution
px = px_bass*0.4+px_salm*0.6;

figure
plot(x,px)
axis([130, 210, 0, inf]);

% posterior probabilities for each fish species
pc_bass = px_bass*0.4./px;
pc_salm = px_salm*0.6./px;

figure

plot(x,pc_bass)
hold on
plot(x,pc_salm)