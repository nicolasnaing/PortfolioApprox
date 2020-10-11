%---------actual numbers----------
TrueCov=[64 0 0; 0 1 0; 0 0 5];
TrueMew=[3; 1.8; 2];

%--------list of possible weights---------
Weights = [];
WeightIndex=1;
for i=0:100
    for j=0:(100-i)
        Weights(1,WeightIndex) = i/100;
        Weights(2,WeightIndex) = j/100;
        Weights(3,WeightIndex) = 1-i/100-j/100;
        WeightIndex=WeightIndex+1;
    end
end

%number to compute is %mewt*x-gamma*xt*cov*x
gam=.2;

%actual
TrueBestVal=0;
TrueBestSeq=[];
for i=1:5151
      val=TrueMew'*Weights(1:3,i)-gam*Weights(1:3,i)'*TrueCov*Weights(1:3,i);
      if val>TrueBestVal
          TrueBestVal=val;
          TrueBestSeq=Weights(1:3,i);
      end
end

TrueBestVal
TrueBestSeq

OtherTrueCov = [0 0 0; 0 0 0; 0 0 0;];

for i=1:3
    for j=1:3
        OtherTrueCov(i,j)=TrueCov(i,j)+TrueCov(j,i);
    end
end

OtherTrueCov;

opt=inv(OtherTrueCov)*TrueMew;
optfinal=opt/gam

opt2=inv(TrueCov)*TrueMew;
optfinal2=opt2/(2*gam)

test = [.12; 9; 1]

return1=TrueMew'*optfinal-gam*optfinal'*TrueCov*optfinal
return2=TrueMew'*optfinal2-gam*optfinal2'*TrueCov*optfinal2
return3=TrueMew'*test-gam*test'*TrueCov*test
%--------generating experimental dataset----- (Est=estimated)
% unif(mean+-2sigma)
% mixture of normals
n=20; %number of time periods for given data
trials=30; %number of trials

for i=1:trials
    EstData=zeros(n,3);
    for i=1:20
        EstData(i,1)=normrnd(3,8);
        EstData(i,2)=normrnd(1.8,1);
        EstData(i,3)=normrnd(2,2)+normrnd(0,1);
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
    EstBestVal=0;
    EstBestSeq=[];
    for i=1:5151
        val=EstMean'*Weights(1:3,i)-gam*Weights(1:3,i)'*EstCovMatrix*Weights(1:3,i);
        if val>EstBestVal
            EstBestVal=val;
            EstBestSeq=Weights(1:3,i);
        end
    end
    EstBestVal;
    EstBestSeq;
    AdjustedValue=TrueMew'*EstBestSeq-gam*EstBestSeq'*TrueCov*EstBestSeq;
    expans=[EstBestVal,AdjustedValue];
    
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
    ShrkBestVal=0;
    ShrkBestSeq=[];
    for i=1:5151
        val=EstMean'*Weights(1:3,i)-gam*Weights(1:3,i)'*ShrunkCov*Weights(1:3,i);
        if val>ShrkBestVal
            ShrkBestVal=val;
            ShrkBestSeq=Weights(1:3,i);
        end
    end
    ShrkBestVal;
    ShrkBestSeq;
    shrkval2=TrueMew'*ShrkBestSeq-gam*ShrkBestSeq'*TrueCov*ShrkBestSeq;
    shrkans=[ShrkBestVal, shrkval2];
end



%trying to use derivatives to find optimal combination; did not go well

OtherMean = [0; 0; 0];
OtherMatrix = [0 0 0; 0 0 0; 0 0 0];

OtherMean(1,1)=EstMean(1)/gam;
OtherMean(2,1)=EstMean(2)/gam;
OtherMean(3,1)=EstMean(3)/gam;


for i=1:3
    for j=1:3
        OtherMatrix(i,j)=EstCovMatrix(i,j)+EstCovMatrix(j,i);
    end
end

OtherMatrix;

inv(OtherMatrix)*OtherMean;

optweights = OtherMean/(OtherMean(1,1)+OtherMean(2,1)+OtherMean(3,1));
optval=EstMean'*optweights-gam*optweights'*ShrunkCov*optweights;

%extra stuff about the pg94 in the book, also did not go well
%lamb=[1 1 1]*inv(EstCovMatrix)*EstMean'/gam;

%firstpart=inv(EstCovMatrix)*EstMean'/gam;
%secondpart = EstCovMatrix*[1; 1; 1;]/([1 1 1]*inv(EstCovMatrix)*[1; 1; 1;]);
%firstpart+(1-lamb)*secondpart




