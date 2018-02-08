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
        if mod(t,round(1e3/(f(n)*T))) ==0
            y1(n,t) = 1;
        end
    end
    trials = [trials n*ones(1,10)];
end

trials = randsample(trials,length(trials));

t=0;
for n = trials
sound(y1(n,:),Tmax/T)
disp('press space')
pause
t = t+1;
answer.E1(t) = input('press 1 if the click sequence sounds discrete, press 0 if it sounds continuous ');

end

save('answerE1.mat','answer')
save('trialsE1.mat','trials')
proba.E1 = zeros(1,length(f));
for n = find(answer.E1 ==1)
    proba.E1(trials(n)) = proba.E1(trials(n))+1;
end

proba.E1 = proba.E1/10;

plot(proba.E1)

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


t=44;
for n = randtrials(45:end)
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
save('randtrials.mat','randtrials')
save('answerE2.mat','answer.E2')
plot(proba.E2(1,:))
hold on 
plot(proba.E2(2,:))
hold on
plot(proba.E2(3,:))





