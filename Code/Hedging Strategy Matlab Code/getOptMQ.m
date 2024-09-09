function [optM,optQ, optLambda, optCV] = getOptMQ(paraSet,C, alpha,sampleAT, sampleZT,sampleDT, K, N)
%implementation of Thm 9
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

paraSetM = paraSet;
paraSetM(3) = 0;
paraSetM(5)= 0;

%%Evaluation lambdaCirc, QCirc and CCirc
lambdaCirc = getLambdaCirc(paraSet, alpha, sampleZT);
disp('lambda circ done.')
[QCirc,CCirc] = getQCCirc(paraSet, alpha, sampleDT, sampleAT, sampleZT, K, N); %[QCirc,CCirc] = getQCCirc(paraSet, al, sampleDT, sampleAT, sampleZT,K,N)
disp('Q circ and C circ done.')

% bmPaths = bmGenerate(N,T);
% [~, sampleDt, sampleAt, sampleZt] = simXtDtAtZt(paraSet,bmPaths,N);
% sampleDT = max(sampleDt(:,end),0);
% sampleAT = sampleAt(:,end);
% sampleZT = sampleZt(:,end);
% sampleDt = [];
% sampleAt = [];
% sampleZt = [];

if C >= CCirc %case (i)
    disp('case i')
    optQ = QCirc;
    optLambda = lambdaCirc; 
    optM = p*optQ - C + 1/mean(sampleZT.*((sampleZT*lambdaCirc)<=1))*max(0, C-CCirc) -K*(optQ > 0);
else  %case (ii)
    disp('case ii')
    [optQ, optM, optLambda] = getCase2Sol(paraSet, C, alpha, QCirc, sampleAT, sampleZT,K);
end



%%calculate opt CV: m - 1/(1-alpha)*sf
%sf = (m-pQ+C)P(optLambda*Z > 1) + (p+c)*E[(hatD^+ \wedge Q - D^+)^+]
sampleHatDT = sampleAT + tilSig*sqrt(T)*norminv(min(1,optLambda*sampleZT));
%sf = (optM - p*optQ + C)*mean((optLambda*sampleZT)>1) + (p+c)*mean(max(0, min(max(sampleHatDT,0),optQ)-max(sampleDT,0)));

%optCV = optM - 1/(1-alpha)*sf ;

cvTerm1 = C*(optLambda/(1-alpha) - 1);
DEvent = (sampleDT < (min(optQ,sampleHatDT))); 
cvTerm2 = (p+c)/(1-alpha)*mean(max(sampleDT,0).*DEvent);

optCV = cvTerm1 + cvTerm2 - K*(optQ>0);

end

%auxilliary function to solve for case (ii)
function [Q2, m2, lambda2] = getCase2Sol(paraSet, C, alpha, QCirc, sampleAT, sampleZT, K)

%implementation of Thm 9
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

paraSetM = paraSet;
paraSetM(3) = 0;
paraSetM(5)= 0;

%%Solve: (p+c)E^M[(Q - \hat D^+_T(lambda))^+] = C and
%        (p+c)E[(lambda * Z) \wedge \Delta(Q)] - p(1-alpha) = 0

%%run bisection on Q
    QL = C/p*0.9;
    QR = QCirc*1.1;
    QM = (QL+QR)/2;
    thisLambda = getLambda2(paraSet, C, QM, sampleAT, sampleZT);
    LHS = (p+c)*mean(min(thisLambda*sampleZT,normcdf((QM-sampleAT)/(tilSig*sqrt(T)))))-p*(1-alpha);
          
    while abs(LHS) > 0.0000000001
        if LHS > 0 %Q too high
            QR = QM;
            QM = (QL+QR)/2;
        else   %Q too low
            QL = QM;
            QM = (QL+QR)/2;
        end
        
        thisLambda = getLambda2(paraSet, C, QM, sampleAT, sampleZT);
        LHS = (p+c)*mean(min(thisLambda*sampleZT,normcdf((QM-sampleAT)/(tilSig*sqrt(T)))))-p*(1-alpha);  
        abs(LHS)
    end
    
    Q2 = QM;
    m2 = p*Q2 - C - K*(Q2>0);
    lambda2 = thisLambda;
    

end

%%auxilliary function to solve for lambda in case (ii)
%(p+c)E^M[(Q - \hat D^+_T(lambda))^+] = C
function thisLambda = getLambda2(paraSet, C,Q, sampleAT, sampleZT)
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

paraSetM = paraSet;
paraSetM(3) = 0;
paraSetM(5)= 0;


%run bisection
lambdaL = 0;
lambdaR = 100;
lambdaM = (lambdaL + lambdaR)/2;
sampleHatD = sampleAT + tilSig*sqrt(T)*norminv(min(1,lambdaM*sampleZT));
LHS = (p+c)*mean(sampleZT.*max(Q - max(sampleHatD,0),0));
iterationCount = 0;
%while abs(LHS-C) > 0.00001 %max(0.000001*C, 0.01)
%while abs(LHS-C)/C > 0.03
while abs(LHS-C) > 0.00001 %max(0.000001*C, 0.01)
    iterationCount = iterationCount + 1; % Increment the iteration counter
    if LHS > C %too low
        lambdaL = lambdaM;
        lambdaM = (lambdaL + lambdaR)/2;
    else %too high
        lambdaR = lambdaM;
        lambdaM = (lambdaL + lambdaR)/2;
    end
    
    sampleHatD = sampleAT + tilSig*sqrt(T)*norminv(min(1,lambdaM*sampleZT));
    LHS = (p+c)*mean(sampleZT.*max(Q - max(sampleHatD,0),0));
    fprintf('Iteration %d: LHS = %.5f, lambdaM = %.5f, C = %.5f, Abs Difference = %.5f\n', iterationCount, LHS, lambdaM, C, abs(LHS - C)); % Print current values to track progress
end
disp('Convergence achieved.');

thisLambda = lambdaM;
end

