function GenerationTimeDistribution=GenerationTimeDistributionFunction

TailThreshold=0.0001;

MeanGamma=3.5;
SDGamma=2.4;

Shape=(MeanGamma/SDGamma)^2;
Scale=MeanGamma/Shape;


GenerationTimeMaxDay=ceil(gaminv(1-TailThreshold,Shape,Scale));
GenerationTimeDistribution=zeros(GenerationTimeMaxDay,1);

i=1;
while i<=GenerationTimeMaxDay
    GenerationTimeDistribution(i)=gamcdf(i,Shape,Scale)-gamcdf(i-1,Shape,Scale);
    i=i+1;
end