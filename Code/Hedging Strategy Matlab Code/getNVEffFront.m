function [nvCVVec, nvMVec, nvQVec, selectionVec, nvKappaVec, nvAlphaVec] = getNVEffFront(paraSet0, paraSet1, alphaVec, sampleD, K, N)
%cvar vs alpha effFront

p0 = paraSet0(1);
c0 = paraSet0(2);
p1 = paraSet1(1);
c1 = paraSet1(2);


T = paraSet1(11);

nvAlphaVec = alphaVec;
nvCVVec = [];
nvMVec = [];
nvQVec = [];
selectionVec = [];
nvKappaVec = [];


for j=1:length(nvAlphaVec)
    disp(['Iteration NV: ', num2str(j)]); % Display current iteration number

    al = alphaVec(j);
    selection = getNVTech(paraSet0, paraSet1,al, K,N, sampleD, T);

    [kappa, ~, ~] = getNVKappa(paraSet0, paraSet1, al, K,N ,T, sampleD);
    nvKappaVec = [nvKappaVec; kappa];
    selectionVec =  [selectionVec;selection];

    if selection == 0
        [Qalpha, nvcv, ~, ~] = getNVCV(paraSet0,al, 0,N, sampleD);
        nvQVec = [nvQVec; Qalpha];
        nvMVec = [nvMVec; p0*Qalpha];
        nvCVVec = [nvCVVec;nvcv];

    else
        [Qalpha, nvcv, ~, ~] = getNVCV(paraSet1,al, K,N, sampleD);
        nvQVec = [nvQVec; Qalpha];
        nvMVec = [nvMVec; p1*Qalpha];
        nvCVVec = [nvCVVec;nvcv];
    end

end