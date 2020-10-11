%---------actual numbers----------
TrueCov=[64 0 0; 0 1 0; 0 0 5];
TrueMew=[3; 1.8; 2];


%number to compute is %mewt*x-gamma*xt*cov*x
gam=.2;
opt2=inv(TrueCov)*TrueMew;
optfinal=opt2/(2*gam)
TrueMew'*optfinal-gam*optfinal'*TrueCov*optfinal

%--------generating experimental dataset----- (Est=estimated)
% unif(mean+-2sigma)
% mixture of normals
n=20; %number of time periods for given data
trials=30; %number of trials

for i=1:170
    EstData=zeros(n,3);
    for i=1:20
        EstData(i,1)=3*rand();
        EstData(i,2)=normrnd(2,2);
        EstData(i,3)=normrnd(3,3);
    end
%----------estimated mean and covariance matrix---------
    EstMean=[0; 0; 0];
    for i=1:n
        EstMean(1,1)=EstMean(1,1)+EstData(i,1)/n;
        EstMean(2,1)=EstMean(2,1)+EstData(i,2)/n;
        EstMean(3,1)=EstMean(3,1)+EstData(i,3)/n;
    end

    EstCovMatrix=[0 0 0; 0 0 0; 0 0 0];

    for i=1:3
        for j=1:3
            xyterm=0;
            xterm=0;
            yterm=0;
            for k=1:n
                xyterm=xyterm+EstData(k,i)*EstData(k,j);
                xterm=xterm+EstData(k,i);
                yterm=yterm+EstData(k,j);
            end
            xyterm=xyterm/n;
            xterm=xterm/n;
            yterm=yterm/n;
            EstCovMatrix(i,j)=xyterm-xterm*yterm;
        end
    end
    
    ExpResult=inv(EstCovMatrix)*EstMean;
    ExpResultFinal=ExpResult/(2*gam);
    FirstReturn=EstMean'*ExpResultFinal-gam*ExpResultFinal'*EstCovMatrix*ExpResultFinal;
    AdjustedReturn=TrueMew'*ExpResultFinal-gam*ExpResultFinal'*TrueCov*ExpResultFinal;
    exppair=[FirstReturn, AdjustedReturn];
    
    
    %experimental with shrunken cov sides
    ShrunkCov=[0 0 0; 0 0 0; 0 0 0];
    for i=1:3
        for j=1:3
            if i==j
                ShrunkCov(i,j)=EstCovMatrix(i,j);
            end
            if i~=j
                ShrunkCov(i,j)=.9*EstCovMatrix(i,j);
            end
        end
    end
    ShrkResult=inv(ShrunkCov)*EstMean;
    ShrkResultFinal=ShrkResult/(2*gam);
    ShrkFirstReturn=EstMean'*ExpResultFinal-gam*ExpResultFinal'*ShrunkCov*ExpResultFinal;
    ShrkAdjustedReturn=TrueMew'*ShrkResultFinal-gam*ShrkResultFinal'*TrueCov*ShrkResultFinal;
    shrkpair9=[ShrkFirstReturn, ShrkAdjustedReturn];
end

