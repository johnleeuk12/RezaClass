function HW4_SBE2()

%% Question 1 

lambda = 100; % spikes per sec
t0 = 2*1e-3; %s


tmax = 0.1; %s
dt = 0.1*1e-3;
timevector = (0:dt:tmax);
p1 = zeros(1,length(timevector));
p2 = zeros(1,length(timevector));
p3 = zeros(1,length(timevector));


for tau = 1:length(timevector)
p1(tau) = lambda*exp(-timevector(tau)*lambda);
end

p2(find(timevector == t0):end-1) = p1(1:end-find(timevector == t0));

F1 = p1/sum(p1);
F2 = p2/sum(p2);

plot(timevector,F1)
hold on
plot(timevector,F2)




%% Question 2 

t02 = 2.5*1e-3;
N1 = 30000;
p3(find(timevector == t02):end-1) = p1(1:end-find(timevector == t02));

C3 = cumsum(p3/sum(p1)); %cdf of poisson process with t0 = t02
ISI1 = [];
for n = 1:N1-1
    tau = randsample(1:length(timevector),1);
    ISI1 = [ISI1 timevector(min(find(ceil(C3*1e3) - tau>=0)))];
    if isempty(find(ceil(C3*1e3) - tau>=0))==1
      ISI1 = [ISI1 timevector(max(find(ceil(C3*1e3) - tau<=0)))];  
    end
end
spikes1 = [0 cumsum(ISI1)];

Count1 = histcounts(ISI1,timevector);
Count1 = Count1/sum(Count1);
plot(Count1)
hold on 
plot(F2)

%% Question 3 

lambda4 = 75;
lambda5 = 55;
N4 = 35000;
N5 = 40000;

p4 = zeros(1,length(timevector));
p5 = zeros(1,length(timevector));
for tau = 1:length(timevector)
p4(tau) = lambda4*exp(-timevector(tau)*lambda4);
p5(tau) = lambda5^2*timevector(tau)*exp(-timevector(tau)*lambda5);
end

C4 = cumsum(p4/sum(p4));
C5 = cumsum(p5/sum(p5));
ISI4 = [];
ISI5 = [];

for n = 1:N4-1
    tau = randsample(1:length(timevector),1);
    ISI4 = [ISI4 timevector(min(find(ceil(C4*1e3) - tau>=0)))];
    if isempty(find(ceil(C4*1e3) - tau>=0))==1
      ISI4 = [ISI4 timevector(max(find(ceil(C4*1e3) - tau<=0)))];  
    end
end

for n = 1:N5-1
    tau = randsample(1:length(timevector),1);
    ISI5 = [ISI5 timevector(min(find(ceil(C5*1e3) - tau>=0)))];
    if isempty(find(ceil(C5*1e3) - tau>=0))==1
      ISI5 = [ISI5 timevector(max(find(ceil(C5*1e3) - tau<=0)))];  
    end
end

spikes3 = sort([0 cumsum(ISI4) cumsum(ISI5)]);
ISI6 = [];
for n = 1:length(spikes3)-1
    ISI6 = [ISI6 spikes3(n+1)-spikes3(n)];
end


Count6 = histcounts(ISI6,timevector);
Count5 = histcounts(ISI5,timevector);

figure
plot(Count6)
hold on 
plot(Count5)

