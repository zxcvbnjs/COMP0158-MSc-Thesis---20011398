function [kappa, qVec, cvVec] = getNVKappa(paraSet0, paraSet1, al, K,N ,T, sampleDT)


p0 = paraSet0(1);
c0 = paraSet0(2);
p1 = paraSet1(1);
c1 = paraSet1(2);

r = p0+c0;

[Q0, cv0, ~, ~] = getNVCV(paraSet0,al,0,N,sampleDT);
[Q1, cv1, ~, ~] = getNVCV(paraSet1,al,K,N,sampleDT);
qVec = [Q0 Q1];
cvVec = [cv0 cv1];


% [~,sampleDT,~] = simDT(paraSet0,T,N);
% sampleDT = max(sampleDT, 0);

kappa = r/(1-al)*mean(sampleDT.*(sampleDT >= Q0).*(sampleDT <= Q1));