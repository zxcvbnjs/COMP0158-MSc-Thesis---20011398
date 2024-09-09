function lambdaCirc = getLambdaCirc(paraSet, alpha, sampleZ)
%E[(1 - \g^\circ Z_T)^+] = \alpha

pp = paraSet(1);
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

% bmPaths = bmGenerate(N,T);
% sampleZt = simZt(paraSet,bmPaths,N);

lambdaL = 0;%(1-alpha)*0.1;
lambdaR = 1000000/alpha;
lambdaM = (lambdaL + lambdaR)/2;
LHS = mean(max(1-lambdaM*sampleZ,0));

disp('Lambda circ bisection begins.')
while abs(LHS - alpha) > 0.00000001
    %LHS-alpha
    
    if LHS > alpha %lambda too low
        lambdaL  = lambdaM;
        lambdaM = (lambdaL + lambdaR)/2;
        LHS = mean(max(1-lambdaM*sampleZ,0));
    else  %lambda too high
        lambdaR = lambdaM;
        lambdaM = (lambdaL + lambdaR)/2;
        LHS = mean(max(1-lambdaM*sampleZ,0));
    end 
end

lambdaCirc = lambdaM;
