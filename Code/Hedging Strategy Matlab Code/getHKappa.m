function [mVec, qVec, lambdaVec, cvVec, kappaH] = getHKappa(paraSet0, paraSet1, alpha, K,N ,T, sampleDT, sampleAT, sampleZT)

p0 = paraSet0(1);
p1 = paraSet1(1);
c0 = paraSet0(2);
c1 = paraSet1(2);

r = p0+c0;

muX = paraSet0(3);
sigX = paraSet0(4);
muE = paraSet0(5);
sigE = paraSet0(6);
rho = paraSet0(7);


tilSig = paraSet0(8);
X0 = paraSet0(9);
E0 = paraSet0(10);

T = paraSet0(11);
C = paraSet0(12);
Z0 = 1;

% Check eta
rhoBar = rho/sqrt(1-rho^2);
etaVec = [muX/sigX rhoBar*(muE/sigE - muX/sigX)];


[optM0,optQ0, optLambda0, optCV0] = getOptMQ(paraSet0,C, alpha,sampleAT, sampleZT,sampleDT, 0, N);
[optM1,optQ1, optLambda1, optCV1] = getOptMQ(paraSet1,C, alpha,sampleAT, sampleZT,sampleDT, K, N);
mVec = [optM0 optM1];
qVec = [optQ0 optQ1];
lambdaVec = [optLambda0 optLambda1];
cvVec = [optCV0 optCV1];


[~,~, sampleDhatT0] = calcVTDhatT(sampleDT,sampleAT,sampleZT,paraSet0,optM0,optQ0,optLambda0,0);
[~,~, sampleDhatT1] = calcVTDhatT(sampleDT,sampleAT,sampleZT,paraSet1,optM1,optQ1,optLambda1,0); %calcVTDhatT(sampleDT,sampleAT,sampleZT,paraSet,m,Q,lambda,gamma)

kappaH = C*(optLambda1 - optLambda0)/(1-alpha)  + r/(1-alpha)*mean(sampleDT.*(sampleDT >= min(optQ0, sampleDhatT0)).*(sampleDT <= min(optQ1, sampleDhatT1)));
