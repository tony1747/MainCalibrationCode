read_CARTdata()
for n_cell=16:18 %check the index
    %for HT1080, the index should be 12
    %for PBT138, the index should be 9
data=cancerdata(:,n_cell);
InteractionCalibration2()
end