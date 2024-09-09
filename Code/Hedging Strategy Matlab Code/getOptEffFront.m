function [optMVec,optQVec, optCVVec, optKappaVec, optSelectionVec, optLambdaVec] = getOptEffFront(paraSet0,paraSet1,alphaVec, sampleA, sampleZ, sampleD, K, N)
kappaH = 100000000;
C = paraSet0(12);
T = paraSet0(11);
optMVec = [];
optQVec = [];
optCVVec = [];
optKappaVec = [];
optSelectionVec = [];
optLambdaVec = [];

for j = 1:length(alphaVec)
    disp(['Iteration CVaR: ', num2str(j)]); % Display current iteration number
    al = alphaVec(j);

    %selection = getOptTech(paraSet0, paraSet1,al, K,N, sampleD,sampleA,sampleZ, T);
    selection = 1;
    %[~, ~, ~, ~, kappaH] = getHKappa(paraSet0, paraSet1, al, K,N ,T, sampleD, sampleA, sampleZ);
    optSelectionVec = [optSelectionVec; selection];
    optKappaVec = [optKappaVec;kappaH];


    if selection ==0 
        [optM,optQ, optLambda, optCV] = getOptMQ(paraSet0,C, al,sampleA, sampleZ,sampleD, 0, N);
        optMVec = [optMVec;optM];
        optQVec = [optQVec; optQ];
        optLambdaVec = [optLambdaVec;optLambda];
        optCVVec = [optCVVec; optCV];
    else
        [optM,optQ, optLambda, optCV] = getOptMQ(paraSet1,C, al,sampleA, sampleZ,sampleD, K, N);
        optMVec = [optMVec;optM];
        optQVec = [optQVec; optQ];
        optLambdaVec = [optLambdaVec;optLambda];
        optCVVec = [optCVVec; optCV];
    end

end




