function RandGenerationTime=RandGenerationTimeFunction(PotentialGenerationTime,PotentialGenerationTimeDistribution)

% rand control function
rng('shuffle');

ProbList=cumsum(PotentialGenerationTimeDistribution);
RandNumber=rand;

RandGenerationTime=PotentialGenerationTime(find(RandNumber<ProbList, 1 ));