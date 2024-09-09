function [bmPaths,sampleDT,A] = simDT(paraSet,T,N)
%paraSet = [p c mu sig tilSig X0]. mu, sig and tilSig are annualized quantities.
%T is the time horizon: 252 points for 1 year, and the unit is 1 year.
%N is the sample size: # DT.
%Note that each DT is computed using one BM and one tilBM


%bmPaths = bmGenerate(N,T); %financial BM
%tilBmPaths = bmGenerate(N,T); %operational BM
%each row is a sample path of length T+1.

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

bmPaths1 = bmGenerate(N,T);
bmPaths2 = bmGenerate(N,T);

bmPathsX = bmPaths1;
bmPathsE = rho*bmPaths1 + sqrt(1-rho^2)*bmPaths2;

bmPaths = [bmPaths1;bmPaths2];
tIncMat = repmat((0:1:(T-0)*252)*(1/252),N,1);

growthMat1 = exp((muX-0.5*sigX^2)*(tIncMat) + sigX*bmPathsX); %the exp(...) term in Xu
growthMat2 = exp((muE-0.5*sigE^2)*(tIncMat) + sigE*bmPathsE); 


rateMat = tilMu(X0*growthMat1, E0*growthMat2);
A = (sum(rateMat(:,2:(end-1)),2)+ (rateMat(:,1)+ rateMat(:,end))/2)*(1/252);
sampleDT = A + tilSig*normrnd(0,sqrt(T),N,1);
   
end