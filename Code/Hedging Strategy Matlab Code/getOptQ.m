function optQ = getOptQ(m,m1, paraSet,N,sampleDT, sampleAT, sampleZT)

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

Q0 = (m+C)/p;

if getQpartialDev(Q0,m,m1,paraSet,N,sampleDT, sampleAT, sampleZT) < 0
    optQ = Q0;
    display('(m + C)/p checked')
elseif getQpartialDev(0,m,m1,paraSet,N,sampleDT, sampleAT, sampleZT) >= 0
    optQ = 0;
   display( '0 checked')
else
    
    % Initilize
    Q_left = 0;
    Q_right = Q0;
    
    numIter = 0;
    
    %tic;
    while abs(Q_right - Q_left)  > 0.1
        
        Q_middle = (Q_left + Q_right)/2;
        y1_middle =  getQpartialDev(Q_middle,m,m1,paraSet,N,sampleDT, sampleAT, sampleZT);
        
        if y1_middle <= 0
            Q_left = Q_middle;
        else
            Q_right = Q_middle;
        end
        
        numIter = numIter + 1;
        %fprintf('Bisection itertation %d\n', numIter);
        %fprintf('Qpartial Derivative = %d\n', y1_middle);
        %fprintf('Error: %d', abs(Q_right - Q_left));
    end
    %toc;
    
    optQ = (Q_left + Q_right)/2;
end

end