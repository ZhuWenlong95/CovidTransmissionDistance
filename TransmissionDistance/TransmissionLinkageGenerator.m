function [IndicatorMatrixSumTheta,IndicatorMatrixSum]=TransmissionLinkageGenerator(OnsetDay)

% rand control function
rng('shuffle');

ThetaNumber=50;

GenerationTimeDistribution=GenerationTimeDistributionFunction;
TransmissionMatrix=TransmissionMatrixFunction(OnsetDay,GenerationTimeDistribution); 
TransmissionMatrix=TransmissionMatrix+TransmissionMatrix';
TransmissionGraph=graph(TransmissionMatrix);

GenerationTimeMaxDay=length(GenerationTimeDistribution);
OnsetDayList=unique(OnsetDay);


IndicatorMatrixSumTheta=zeros(OnsetDayList(end),OnsetDayList(end),ThetaNumber);
IndicatorMatrixSum=zeros(OnsetDayList(end),OnsetDayList(end));

i=1;
while i<=OnsetDayList(end)-1
    j=i+1;
    while j<=i+GenerationTimeMaxDay && j<=OnsetDayList(end)
        IndicatorMatrixSum(i,j)=sum(OnsetDay==i)*sum(OnsetDay==j);
        List1=find(OnsetDay==i);
        List2=find(OnsetDay==j);
        k=1;
        while k<=length(List1)
            h=1;
            while h<=length(List2)
                ShortestDistance=distances(TransmissionGraph,List1(k),List2(h));
                if isfinite(ShortestDistance) && ShortestDistance<=ThetaNumber
                    IndicatorMatrixSumTheta(i,j,ShortestDistance)=IndicatorMatrixSumTheta(i,j,ShortestDistance)+1;
                end
                h=h+1;
            end
            k=k+1;
        end
        j=j+1;
    end
    i=i+1;
end