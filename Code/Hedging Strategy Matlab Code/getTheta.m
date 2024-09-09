function thetaVec = getTheta(paraSet,m,Q,lambda,t,X,A,N)

%Refer to the paper for theta.

%X and A are values of Xt and At at time t.

p = paraSet(1);
c = paraSet(2);
mu = paraSet(3);
sig = paraSet(4);
tilSig = paraSet(5);
X0 = paraSet(6);
T = paraSet(7); 
C = paraSet(8);
Z0 = 1;

eta = mu/sig;

% tau(u) = u - t;
tIncVec = [0:1:((T-t)*252)]*(1/252);
tIncMat = repmat(tIncVec,N,1);
bmIncMat = bmIncGenerate(t,T,N);

K = (m + C + c*Q)/(p+c);

Z = exp(-0.5*(mu - eta^2)*t)*(X/X0)^(-eta^2/mu);
%Z

xCircMat = exp(-0.5*sig^2*tIncMat + sig*bmIncMat);
zCirc = exp(0.5*eta^2*(T-t) - eta*bmIncMat(:,end));

%rateMat = tilMu(Xt(k)*growthMat);
%intVec = (sum(rateMat(:,2:(end-1)),2)+(rateMat(:,1)+ rateMat(:,end))/2)*(1/252); %trapezoidal rule

integrandMat2 = tilMuDev(X*xCircMat).*xCircMat;
integrandMat1 = tilMu(X*xCircMat);

%integralPart2 is the integral involving \tilde\mu^\prime
%integralPart1 is the integral to construct A_T.
integralPart2 = (sum(integrandMat2(:,2:(end-1)),2) + (integrandMat2(:,1) + integrandMat2(:,end))/2)*(1/252);
integralPart1 = (sum(integrandMat1(:,2:(end-1)),2) + (integrandMat1(:,1) + integrandMat1(:,end))/2)*(1/252);
%mean(integralPart2)


%define theta1 and theta2 below.
theta1 = (p+c)*max(K - Q, 0)/(sig*X*sqrt(T - t))*normpdf(log(lambda*Z)/(eta*sqrt(T-t)) - 0.5*eta*sqrt(T - t));
%normpdf(log(lambda*Z)/(eta*sqrt(T-t)) - 0.5*eta*sqrt(T - t))
%for theta2, need to use \hat D and Z_T.
zVec = Z*zCirc;
hatDVec = A + integralPart1 + tilSig*sqrt(T)*norminv(min(lambda*zVec,1));

%mean(hatDVec)

%indicator vec
theta2Vec1 = -(p + c)*((hatDVec >=0).*(hatDVec <= min(K,Q)));
theta2Vec2 = integralPart2 - tilSig*sqrt(T)*(eta/(sig*X))*lambda*zVec./normpdf(norminv(min(1,lambda*zVec)));
%Re1(isinf(Re1)) = 0;
%theta2Vec2(isinf(theta2Vec2)) = 0;
theta2Vec2(hatDVec > min(K,Q)) = 0;

theta2 = mean(theta2Vec1.*theta2Vec2);

%mean(normpdf(norminv(min(1,lambda*zVec))))
%mean(tilSig*sqrt(T)*(eta/(sig*X))*lambda*zVec./normpdf(norminv(min(1,lambda*zVec))))

thetaVec = [theta1;theta2];
end