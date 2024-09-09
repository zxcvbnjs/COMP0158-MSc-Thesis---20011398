function [Qalpha, nvcv, meanH, varH] = getNVCV(paraSet,al, K,N, sampleDT)
%paraSet = [p c mu sig tilSig X0]. mu, sig and tilSig are annualized quantities.
%this function computes the Q^* in the calssical NV model (hence in the physical measure), and the
%corresponding optimal operational payoff


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

beta = p/(p+c)*(1-al);

%bmPaths = bmGenerate(N,T);
%[~, sampleDt] = simXtDt(paraSet,bmPaths,T,N); %A set of sample for DT is generated.

%[~,sampleDT,~] = simDT(paraSet,T,N);
%sampleDT = max(sampleDT, 0);

%sampleDT = sampleDt(:,end);
%display('Simulation done')
%Now get the inverse cdf function by bisection method:
lb = 10000;
ub = 10000000;
midpoint = (lb+ub)/2;
 currentLoc = getEmDCDF(sampleDT,midpoint);
while  abs(currentLoc-beta) > 1/N
    if currentLoc < beta %lower bound should be updated
        lb = midpoint;
    else %upper bound should be updated
        ub = midpoint;
    end
    midpoint = (lb+ub)/2;
    currentLoc = getEmDCDF(sampleDT,midpoint);
    fprintf('The currentLoc is: %d\n', currentLoc);
end

Qalpha = midpoint;
hSample = p*Qalpha - (p+c)*max(Qalpha-sampleDT,0) - K*(Qalpha > 0);
meanH = sum(hSample)/length(hSample);
varH = var(hSample);

%optM = p*Qalpha - K*(Qalpha > 0);
%U = m - \frac{1}{1-\alpha}\ex[(m - H)^+]
sampleU = (p+c)/(1-al)*mean(sampleDT.*(sampleDT <= Qalpha)) - K*(Qalpha > 0);
nvcv = mean(sampleU);
end


function emHCDF = getEmDCDF(sampleDT,x)
%this computes the empiracle cdf value at x, for a *given* sample
N = length(sampleDT);

emHCDF = sum(sampleDT<=x)/N;
end




