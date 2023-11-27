clc
clear

load MeanDistanceShanghai
load RawDataShanghai


Scale=5;
%200M -> 0

MeanDistance=MeanDistance*Scale;
Lamda=1/(MeanDistance-1)+2;
DistanceList=(1:15)'*Scale;
% https://wiki.analytica.com/index.php?title=Power_law_distribution#DensPowerLawDist.28_x.2C_lambda.29
plist=DistanceList.^(1-Lamda);
