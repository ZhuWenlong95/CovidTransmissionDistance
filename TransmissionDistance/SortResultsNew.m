

clc
clear


ResultsFolder='results';
MeanTransmissionFolder='MeanTransmissionFolder';
if ~exist(MeanTransmissionFolder,'dir')
    mkdir(MeanTransmissionFolder)
end

cd RawData
RawDataFile=dir('*.mat');
cd ..

h=1;
while h<=length(RawDataFile)
    
    cd RawData
    load(RawDataFile(h).name)  
    cd ..      
    cd(ResultsFolder)
    load(['IndicatorMatrix',RawDataFile(h).name])  
    cd ..

    ThetaNumber=50;

    OmegaMatrix=zeros(size(IndicatorMatrixSumTotal,1),size(IndicatorMatrixSumTotal,2),ThetaNumber);

    TempSum=sum(IndicatorMatrixSumThetaTotal,4);
    i=1;
    while i<=ThetaNumber
        OmegaMatrix(:,:,i)=TempSum(:,:,i)./IndicatorMatrixSumTotal(:,:,1)/size(IndicatorMatrixSumThetaTotal,4);
        i=i+1;
    end
    OmegaMatrix(isnan(OmegaMatrix))=0;

    OnsetDayList=unique(OnsetDay);
    MeanDistanceMatrix=zeros(size(IndicatorMatrixSumTotal,1),size(IndicatorMatrixSumTotal,2));
    i=1;
    while i<=OnsetDayList(end)
        j=i+1;
        while j<=OnsetDayList(end)
            TempListi=find(OnsetDay==i);
            TempListj=find(OnsetDay==j);
            TempDistance=0;
            k=1;
            while k<=length(TempListi)
                l=1;
                while l<=length(TempListj)
                    TempDistance=TempDistance+DistanceEquation([LatitudeData(TempListi(k)),LongitudeData(TempListi(k))],[LatitudeData(TempListj(l)),LongitudeData(TempListj(l))]);
                    l=l+1;
                end
                k=k+1;
            end
            if length(TempListi)*length(TempListj)>0
                MeanDistanceMatrix(i,j)=TempDistance/length(TempListi)/length(TempListj);
            end
            j=j+1;
        end
        i=i+1;
    end


    MeanTransmissionMatrix=zeros(size(IndicatorMatrixSumTotal,1),size(IndicatorMatrixSumTotal,2));

    i=1;
    while i<=OnsetDayList(end)
        j=i+1;
        while j<=OnsetDayList(end)
            TempOmegaSum=0;
            k=1;
            while k<=ThetaNumber
                TempOmegaSum=TempOmegaSum+OmegaMatrix(i,j,k)*sqrt(2*pi*k);
                k=k+1;
            end
            if TempOmegaSum>0
                MeanTransmissionMatrix(i,j)=2*MeanDistanceMatrix(i,j)/TempOmegaSum;
            end
            j=j+1;
        end
        i=i+1;
    end

    MeanDistanceSum=0;
    LinkSum=0;
    i=1;
    while i<=OnsetDayList(end)
        j=i+1;
        while j<=OnsetDayList(end)
            TempListi=find(OnsetDay==i);
            TempListj=find(OnsetDay==j);
            MeanDistanceSum=MeanDistanceSum+MeanTransmissionMatrix(i,j)*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        i=i+1;
    end
    MeanDistance=MeanDistanceSum/LinkSum;


    MeanDistanceVarianceSum=0;
    LinkSum=0;
    i=1;
    while i<=OnsetDayList(end)
        j=i+1;
        while j<=OnsetDayList(end)
            TempListi=find(OnsetDay==i);
            TempListj=find(OnsetDay==j);
            MeanDistanceVarianceSum=MeanDistanceVarianceSum+(MeanTransmissionMatrix(i,j)-MeanDistance)^2*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        i=i+1;
    end
    MeanDistanceVariance=MeanDistanceVarianceSum/(LinkSum-1);
    ConnectionLinkNumber=LinkSum;


    TempList=sum(MeanDistanceMatrix>0,2);
    TempList=find(TempList>1);

    MeanDistanceDate=zeros(TempList(end),1);
    MeanDistanceDateVariance=zeros(TempList(end),1);
    i=1;
    while i<=length(MeanDistanceDate)
        MeanDistanceSum=0;
        LinkSum=0;
        j=i+1;
        while j<=size(MeanDistanceMatrix,1)
            TempListi=find(OnsetDay==i);
            TempListj=find(OnsetDay==j);
            MeanDistanceSum=MeanDistanceSum+MeanTransmissionMatrix(i,j)*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        MeanDistanceDate(i)=MeanDistanceSum/LinkSum;
        MeanDistanceSum=0;
        LinkSum=0;
        j=i+1;
        while j<=size(MeanDistanceMatrix,1)
            TempListi=find(OnsetDay==i);
            TempListj=find(OnsetDay==j);
            MeanDistanceSum=MeanDistanceSum+(MeanTransmissionMatrix(i,j)-MeanDistanceDate(i))^2*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        MeanDistanceDateVariance(i)=MeanDistanceSum/(LinkSum-1);
        i=i+1;
    end


    MeanDistanceDate-norminv(1-0.05/2)*sqrt(MeanDistanceDateVariance)
    MeanDistanceDate+norminv(1-0.05/2)*sqrt(MeanDistanceDateVariance)

    MeanDistanceTime=zeros(14,1);
    i=1;
    while i<=length(MeanDistanceTime)
        MeanDistanceSum=0;
        LinkSum=0;
        j=1;
        while i+j<=size(MeanDistanceMatrix,1)
            TempListi=find(OnsetDay==j);
            TempListj=find(OnsetDay==j+i);
            MeanDistanceSum=MeanDistanceSum+MeanTransmissionMatrix(j,j+i)*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        MeanDistanceTime(i)=MeanDistanceSum/LinkSum;
        i=i+1;
    end


    TempList=sum(MeanDistanceMatrix>0,2);
    TempList=find(TempList>1);
    MeanDistanceDateEND=zeros(TempList(end),1);
    MeanDistanceDateVariance=zeros(TempList(end),1);
    i=1;
    while i<=length(MeanDistanceDateEND)
        MeanDistanceSum=0;
        LinkSum=0;
        j=1;
        while j<=i
            TempListi=find(OnsetDay==i+1);
            TempListj=find(OnsetDay==j);
            MeanDistanceSum=MeanDistanceSum+MeanTransmissionMatrix(j,i+1)*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        MeanDistanceDateEND(i)=MeanDistanceSum/LinkSum;
        MeanDistanceSum=0;
        LinkSum=0;
        j=1;
        while j<=i
            TempListi=find(OnsetDay==i+1);
            TempListj=find(OnsetDay==j);
            MeanDistanceSum=MeanDistanceSum+(MeanTransmissionMatrix(j,i+1)-MeanDistanceDateEND(i))^2*length(TempListi)*length(TempListj);
            LinkSum=LinkSum+length(TempListi)*length(TempListj);
            j=j+1;
        end
        MeanDistanceDateVariance(i)=MeanDistanceSum/(LinkSum-1);
        i=i+1;
    end
    
    cd(MeanTransmissionFolder)
        TempFileName=['MeanTransmission',RawDataFile(h).name];
        save(TempFileName,'MeanTransmissionMatrix','MeanDistance','MeanDistanceTime','MeanDistanceDate','MeanDistanceVariance','ConnectionLinkNumber','MeanDistanceDateEND','-v6');
    cd ..
    
    h=h+1;
end
    
% MeanDistance-norminv(1-0.05/2)*sqrt(MeanDistanceVariance)
% MeanDistance+norminv(1-0.05/2)*sqrt(MeanDistanceVariance)
% 
% 
% MeanDistanceDateEND-norminv(1-0.05/2)*sqrt(MeanDistanceDateVariance)
% MeanDistanceDateEND+norminv(1-0.05/2)*sqrt(MeanDistanceDateVariance)