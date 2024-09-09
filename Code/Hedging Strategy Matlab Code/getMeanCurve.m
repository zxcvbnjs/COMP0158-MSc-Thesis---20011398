function meanCurve = getMeanCurve(paraSet,N,qCurve,QStar,lambdaCurve)

%qCurve = [mVec QVec]
%lambdaCurve = [mVec lambdaVec]
%meanCurve = [mVec meanW meanV meanNVBest wLB]
%wLB is the LB on E[W] in (99) of v5c

p = paraSet(1);
c = paraSet(2); 
mu = paraSet(3);
sig = paraSet(4);
tilSig = paraSet(5);
X0 = paraSet(6);
T = paraSet(7);
C = paraSet(8);
Z0 = 1;

meanCurve = zeros(size(qCurve,1),5);
meanCurve(:,1) = qCurve(:,1);

for k = 1:size(qCurve,1)
    k
bmPaths = bmGenerate(N,T);
[~, sampleDt, sampleAt, sampleZt] = simXtDtAtZt(paraSet,bmPaths,N);

%[WT,VT, DhatT] = calcVTDhatT(sampleDt,sampleAt,sampleZt,paraSet,m,Q,lambda,gamma)
[WT,VT, ~] = calcVTDhatT(sampleDt,sampleAt,sampleZt,paraSet,qCurve(k,1),qCurve(k,2),lambdaCurve(k,2),0);
%HT = p*qCurve(k,2) - (p + c)*max(qCurve(k,2) - max(sampleDt(:,end),0),0);

%Q = (m/p) \wedge Q^*
nvQ = min(qCurve(k,1)/p,QStar);
nvHT = p*nvQ - (p + c)*max(nvQ - max(sampleDt(:,end),0),0);

meanCurve(k,2) = mean(WT);
meanCurve(k,3) = mean(VT);
meanCurve(k,4) = mean(nvHT);

%Get the LB
[~,psi] = getBetaPsi(paraSet,lambdaCurve(k,2));
meanCurve(k,5) = meanCurve(k,4) + C*psi;


end


end