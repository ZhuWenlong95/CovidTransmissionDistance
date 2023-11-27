clc
clear

TrialTimes=500;

load GAMMAParameter

cd RawData
delete ._*
RawDataFile=dir('*.mat');
cd ..

i=1;
while i<=length(RawDataFile)
    cd RawData
    load(RawDataFile(i).name)
    cd ..
    
    
    ThetaNumber=50;

    while sum(ConfirmedList>(1:length(ConfirmedList))')>0
        TempList=find(ConfirmedList>(1:length(ConfirmedList))');
        TempNumber1=TempList(1);
        TempNumber2=ConfirmedList(TempNumber1);
        LatitudeData=LatitudeData([1:TempNumber1-1,TempNumber2,TempNumber1+1:TempNumber2-1,TempNumber1,TempNumber2+1:end]);
        LongitudeData=LongitudeData([1:TempNumber1-1,TempNumber2,TempNumber1+1:TempNumber2-1,TempNumber1,TempNumber2+1:end]);
        OnsetDay=OnsetDay([1:TempNumber1-1,TempNumber2,TempNumber1+1:TempNumber2-1,TempNumber1,TempNumber2+1:end]);
        ConfirmedList=ConfirmedList([1:TempNumber1-1,TempNumber2,TempNumber1+1:TempNumber2-1,TempNumber1,TempNumber2+1:end]);
        ConfirmedList(ConfirmedList==TempNumber2)=TempNumber1;
    end
    
    

    OnsetDayList=unique(OnsetDay);

    IndicatorMatrixSumThetaTotal=zeros(OnsetDayList(end),OnsetDayList(end),ThetaNumber,TrialTimes);
    IndicatorMatrixSumTotal=zeros(OnsetDayList(end),OnsetDayList(end),TrialTimes);
    MeanGamma=GAMMAParameter(str2double(RawDataFile(i).name(1)),1);
    SDGamma=GAMMAParameter(str2double(RawDataFile(i).name(1)),2);

    delete(gcp('nocreate'))
    p=gcp;
    parpoolNumber=round(p.NumWorkers);
    delete(gcp('nocreate'))
    parpool(parpoolNumber);
    parfor h=1:TrialTimes
        [IndicatorMatrixSumTheta,IndicatorMatrixSum]=TransmissionLinkageGenerator3(OnsetDay,ConfirmedList,MeanGamma,SDGamma);
        IndicatorMatrixSumThetaTotal(:,:,:,h)=IndicatorMatrixSumTheta;
        IndicatorMatrixSumTotal(:,:,h)=IndicatorMatrixSum;
    end
    FileName=['IndicatorMatrix',RawDataFile(i).name];
    save(FileName,'IndicatorMatrixSumThetaTotal','IndicatorMatrixSumTotal');
    
    i=i+1;
end

