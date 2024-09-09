function nvSFCurve = getNVSFCurve(paraSet,mVec,QStar,N)

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

sfVec = [];
bmPaths1 = bmGenerate(N,T);
bmPaths2 = bmGenerate(N,T);
[~, ~, sampleDt, ~, ~] = simXtEtDtAtZt(paraSet,bmPaths1,bmPaths2, N);

for k = 1:length(mVec)
    thisQNV = min(mVec(k)/p, QStar);
    
    %Big bug found on April 19 2018; not affected existing results;
    %Because for m \leq p\nvq, Q = m/p, rendering [m + cQ - (p+c)D^+]^+ 
    %= m-pQ + (p+c)(Q-D+)^+ = (p+c)(Q-D^+)^+
    %thisSFVec = max(mVec(k) + c*thisQNV - (p + c)*max(sampleDt(:,end),0),0);
    
    %corrected
    thisSFVec = max(mVec(k) - p*thisQNV+(p+c)*max(thisQNV - max(sampleDt(:,end),0),0),0);
    
    sfVec = [sfVec;mean(thisSFVec)];    
end

nvSFCurve = [mVec sfVec];


end