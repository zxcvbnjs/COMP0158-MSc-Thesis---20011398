function lambdaQCurve = getLambdaQCurve(paraSet,m,m1,qVec,N)
%For a given m, plot \g(m,Q) against Q.


lambdaQCurve = zeros(length(qVec),2);
lambdaQCurve(,:) = qVec;


for k = 1:length(qVec)
    
    lambdaQCurve(k,2) = getLambda(paraSet,m,m1,qVec(k),N);
    
end


