TrueCov=[.75 0 0; 0 4 0; 0 0 9];
TrueMew=[1.5; 2; 3];


%number to compute is %mewt*x-gamma*xt*cov*x
gam=.2;
opt2=inv(TrueCov)*TrueMew;
optfinal=opt2/(2*gam);
TrueMew'*optfinal-gam*optfinal'*TrueCov*optfinal;

results=zeros(200,20);

%--------generating experimental dataset----- (Est=estimated)
% unif(mean+-2sigma)
% mixture of normals
n=20; %number of time periods for given data
trials=30; %number of trials

for f=1:1:10
   g=-.01+f/100;
   for z=1:200
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
        %experimental with shrunken cov sides
        ShrunkCov=[0 0 0; 0 0 0; 0 0 0];
        for i=1:3
            for j=1:3
                if i==j
                    ShrunkCov(i,j)=EstCovMatrix(i,j);
                end
                if i~=j
                    ShrunkCov(i,j)=g*EstCovMatrix(i,j);
                end
            end
        end
        ShrkResult=inv(ShrunkCov)*EstMean;
        ShrkResultFinal=ShrkResult/(2*gam);
        ShrkFirstReturn=EstMean'*ShrkResultFinal-gam*ShrkResultFinal'*ShrunkCov*ShrkResultFinal;
        ShrkAdjustedReturn=TrueMew'*ShrkResultFinal-gam*ShrkResultFinal'*TrueCov*ShrkResultFinal;
        results(z,2*f-1)=ShrkFirstReturn;
        results(z,2*f)=ShrkAdjustedReturn;
   end
end
results