function [mcirc,undm, ovrm] = getMCirc(paraSet,N)

%Implement Lemma 9(ii) of v5c.

p = paraSet(1);
c = paraSet(2);
mu = paraSet(3);
sig = paraSet(4);
tilSig = paraSet(5);
X0 = paraSet(6);
T = paraSet(7);
C = paraSet(8);
Z0 = 1;

% Check eta
eta = mu/sig;


undm = getUndm(paraSet);
[ovrm,~] = getOvrm(paraSet,N);

mcirc = max(undm,ovrm);