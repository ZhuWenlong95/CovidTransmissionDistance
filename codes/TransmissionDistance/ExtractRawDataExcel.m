clc
clear

cd RawDataExcel
delete ._*
DataList=dir('*.xlsx');
cd ..


i=1;
while i<=length(DataList)
    cd RawDataExcel
    RawTable = readtable(DataList(i).name);
    cd ..
    cd RawData
    ConfirmedList=zeros(size(RawTable,1),1);
    LatitudeData=table2array(RawTable(:,"lat"));
    LongitudeData=table2array(RawTable(:,"lng"));
    Temp=datenum((table2array(RawTable(:,"date"))));
    OnsetDay=Temp-Temp(1)+1;
    if length(RawTable.Properties.VariableNames)>4
        Relation=table2array(RawTable(:,"relation"));
        ID=table2array(RawTable(:,"id"));
        j=1;
        while j<=length(Relation)
            Temp=find(strcmpi(ID,Relation{j}));
            if ~isempty(Temp)
                ConfirmedList(j)=Temp;
            end            
            j=j+1;
        end
    end
    save([DataList(i).name(1:end-5),'.mat'],'ConfirmedList','LatitudeData','LongitudeData','OnsetDay');
    cd ..
    i=i+1;
end