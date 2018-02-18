function HW1SB2()

clear all
%% Experiment 1

f = 10:10:100; %repetition rates, in Hz
y1 = zeros(length(f),1e4); % y is a series of vectors used as input to create click trains
% ICI = 1/f*1e3; %ms

Tmax = 1e3; %ms
T = 0.1; %ms
trials = [];
for n = 1:length(f) 
    for t = 1:Tmax/T
        if mod(t,round(1e3/(f(n)*T))) ==0 % creating clicktrains for different repetition rates. 
            y1(n,t) = 1;
        end
    end
    trials = [trials n*ones(1,10)]; %creating a vector for trial number randomization. 
end

trials = randsample(trials,length(trials));

t=0;
for n = trials %playing sounds in a randomized way. 
sound(y1(n,:),Tmax/T)
disp('press space')
pause
t = t+1;
answer.E1(t) = input('press 1 if the click sequence sounds discrete, press 0 if it sounds continuous ');

end

% % save('answerE1.mat','answer')
% % save('trialsE1.mat','trials')
proba.E1 = zeros(1,length(f)); %calculating probability of responding '1'.
for n = find(answer.E1 ==1)
    proba.E1(trials(n)) = proba.E1(trials(n))+1;
end

proba.E1 = proba.E1/10;



figure
plot(f,proba.E1)
title('psychometric function for flutter/fusion discrimination')
ylabel('probability')
xlabel('Stimulus Frequency (Hz)')

%% Experiment 2

ICI = [10 30 80];%ms
jitter = [5:5:25]*1e-2;
nbtrials = 20;
trials = [];
for j = 1:length(jitter)
    trials = [trials ones(1,nbtrials)*j];
end
for ipi = 1:length(ICI)
    y2{ipi} = zeros(length(jitter)*nbtrials,1e4);
    for tr = 1:length(trials)
        for t = 1:Tmax/T
            if mod(t,round(ICI(ipi)/T))==0
                
                tmod = t+ round(ICI(ipi)/T*rand*jitter(trials(tr))*randsample([-1 1],1)/T);
                if tmod <Tmax/T && tmod>1
                y2{ipi}(tr,tmod) = 1;
                end
            end
        end
    end
end

randtrials = randsample(trials,length(trials));
answer.E2=zeros(3,length(randtrials));


t=0
for n = randtrials
    disp('press space')
    pause
    t = t+1
    sound(y2{3}(n,:),Tmax/T)
    
    answer.E2(3,t) = input('press 1 if the click sequence sounds aperiodic, press 0 if it sounds periodic ');
end

proba.E2 = zeros(3,length(jitter));
for ipi = 1:length(ICI)
    for n = find(answer.E2(ipi,:) ==1)
        proba.E2(ipi,randtrials(n)) = proba.E2(ipi,randtrials(n))+1;
    end
end


proba.E2 = proba.E2/20;
% save('randtrials.mat','randtrials')
figure
plot(jitter,proba.E2(1,:))
hold on 
plot(jitter,proba.E2(2,:))
hold on
plot(jitter,proba.E2(3,:))

title('psychometric  functions for periodicity')

xlabel('jitter')
ylabel('probability of the stimulus being perceived as aperiodic')
legend('10ms','30ms', '80ms')



%% Experiment 3

X = 1:0.5:20;

clickA = unifpdf(X,7.5,12.5);
clickB = unifpdf(X,10,14);
clickC = unifpdf(X,11,15);
figure
plot(X,clickA)
hold on 
plot(X,clickB)
legend('click train A','click train B')
title('ICI distribution of two click trains')
xlabel('ICI distribution (ms)')
axis([0 20 0 0.5])

%Calculate probabilities of hits and false alarms


beta = 6:15; %ms
Pd = zeros(1,length(beta));
Pd2 = zeros(1,length(beta));
Pfa = zeros(1,length(beta));

for b = 1:length(beta)
    if beta(b)<= 10
        Pd(b) = 1;
    elseif beta(b)>=14
        Pd(b) = 0;
    else
        Pd(b) = (14-beta(b))*clickB(24);
    end
    if beta(b)<= 11
        Pd2(b) = 1;
    elseif beta(b)>=15
        Pd2(b) = 0;
    else
        Pd2(b) = (15-beta(b))*clickC(24);
    end
    if beta(b)<=7.5
        Pfa(b) = 1;
    elseif beta(b)>=12.5
        Pfa(b) = 0;
    else
        Pfa(b) = (12.5-beta(b))*clickA(20);
    end
end


figure
plot(Pfa,Pd)
hold on
plot(Pfa,Pd2)
title('ROC curve')
ylabel('Probability of hits')
xlabel('probability of false alarms')
legend('Click train B [10 14]','Click train B [11 15]');













