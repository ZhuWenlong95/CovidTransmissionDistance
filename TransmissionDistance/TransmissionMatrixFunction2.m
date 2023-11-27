function TransmissionMatrix=TransmissionMatrixFunction2(OnsetDay,GenerationTimeDistribution,ConfirmedList) 

% rand control function
rng('shuffle');

SparseI=[];
SparseJ=[];
SparseK=[];

GenerationTimeMaxDay=length(GenerationTimeDistribution);
OnsetDayList=unique(OnsetDay);
i=length(OnsetDayList);
while i>1
    StudyDay=OnsetDayList(i);
    PotentialGenerationTime=sort(StudyDay-intersect(OnsetDayList,StudyDay-1:-1:StudyDay-GenerationTimeMaxDay));
    PotentialGenerationTimeDistribution=GenerationTimeDistribution(PotentialGenerationTime)/sum(GenerationTimeDistribution(PotentialGenerationTime));
    StudyList=find(OnsetDay==StudyDay);
    j=1;
    while j<=length(StudyList)
        PotentialGenerationTimePerson=RandGenerationTimeFunction(PotentialGenerationTime,PotentialGenerationTimeDistribution);
        SparseJ=[SparseJ,StudyList(j)];
        SparseK=[SparseK,1];
        InfectList=find(OnsetDay==StudyDay-PotentialGenerationTimePerson);
        SparseI=[SparseI,InfectList(ceil(rand*length(InfectList)))];
        j=j+1;
    end
    i=i-1;
end

TempList=find(ConfirmedList>0);
i=1;
while i<=length(TempList)
    TempPosi=find(SparseJ==TempList(i));
    SparseI(TempPosi)=ConfirmedList(TempList(i));
    i=i+1;
end

TransmissionMatrix = sparse(SparseI,SparseJ,SparseK,length(OnsetDay),length(OnsetDay));