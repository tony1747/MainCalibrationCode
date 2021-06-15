

for receptor = ['L', 'M', 'H', 'VH']
  for CARTnum = [5, 10, 20]
    read_CARTdata() %% <- change the last part, data to cancerdata <- maybe you can also make something to select the cancer cell type outside. 
    CARTratio=1/CARTnum;
    for n = 1:3
      data = selected_cancerdata(:,n); 
      InteractionCalibration_slow(); 
    end
    
    %% here save parameters/posterior 
  end
end 