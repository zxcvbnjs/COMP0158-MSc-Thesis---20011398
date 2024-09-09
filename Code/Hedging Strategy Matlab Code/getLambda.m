function lambda = getLambda(paraSet,m,m1,Q,N, sampleDT, sampleAT, sampleZT)
% This is for short fall model *without* mean constraint

% Stopping Criteria: lambda_right - lambda_left < 0.0000001

gamma = 0;

% Initilize
lambda_left = 0;
lambda_right = 10;

numIter = 0;
%tic;
while abs(lambda_right - lambda_left) > 10^-8
     
   %abs(lambda_right - lambda_left)
    lambda_middle = (lambda_left + lambda_right)/2;
    [y1_middle, ~] = KKTfunc([lambda_middle gamma], paraSet,m,m1,Q,N, sampleDT, sampleAT, sampleZT);
    
    if y1_middle > 0
        lambda_left = lambda_middle;
    else
        lambda_right = lambda_middle;
    end
    
    numIter = numIter + 1;
    %fprintf('Itertation %d\n', numIter);
    %fprintf('y1 = %d\n', y1_middle);
end
%toc;

lambda = (lambda_left + lambda_right)/2;
end