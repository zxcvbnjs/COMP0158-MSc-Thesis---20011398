function [y1, y2] = KKTfunc(multiplierVec, paraSet,m,m1,Q,N, sampleDT, sampleAT, sampleZT)

%m is SF-target
%m1 is mean-constraint
%multiplierVec = [lambda; gamma]




lambda = multiplierVec(1);
gamma = multiplierVec(2);



p = paraSet(1);
c = paraSet(2);
muX = paraSet(3);
sigX = paraSet(4);
muE = paraSet(5);
sigE = paraSet(6);
rho = paraSet(7);


tilSig = paraSet(8);
X0 = paraSet(9);
E0 = paraSet(10);

T = paraSet(11);
C = paraSet(12);
Z0 = 1;

% Check eta
rhoBar = rho/sqrt(1-rho^2);
etaVec = [muX/sigX rhoBar*(muE/sigE - muX/sigX)];

%a = (-0.5*eta^2*T - log(lambda))/(eta*sqrt(T));
a = -(log(lambda) + 0.5*(norm(etaVec))^2*T)/(norm(etaVec)*sqrt(T));


%bmPaths1 = bmGenerate(N,T);
%bmPaths2 = bmGenerate(N,T);


%[~, ~, sampleDt, sampleAt, sampleZt] = simXtEtDtAtZt(paraSet,bmPaths1,bmPaths2, N); 
exH = p*Q - (p + c)*mean(max(Q - max(sampleDT,0),0));
delta = m1 - exH;
%[WT,VT, DhatT] = calcVTDhatT(sampleDt,sampleAt,sampleZt,paraSet,m,Q,lambda,gamma)
[~,VT, DhatT] = calcVTDhatT(sampleDT,sampleAT,sampleZT,paraSet,m,Q,lambda,gamma);

%max((Q - max(DhatT,0)),0))
putOpt = (p + c)*max(Q - max(DhatT,0),0);

ZT = sampleZT;

if lambda > 0 && gamma >= 0
    
    y1 = (m - p*Q) - (m - p*Q + C)*normcdf(a,'upper') + mean(putOpt.*ZT);  %mean(ZT .* VT);
    %y1 = (m - p*Q + C)*normcdf(a)/C +  mean(putOpt.*ZT)/C - 1;
    y2 = gamma * (mean(VT) - delta);
else
    y1 = 100;
    y2 = 100;
end

end