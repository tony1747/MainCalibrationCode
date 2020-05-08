%% code to read excel files 

cancerdata = xlsread( 'Run1_HT1080.xlsx', 'Sheet1' ); 
%cancerdata = xlsread( 'Run8_PBT138.xlsx', 'Sheet1' ); 


% Time in first column (in sec) 
time = cancerdata(:,1); 
time = time/3600 /24;
ind = find(time>1.25,1);
time = time(ind:end);

cancerdata  = cancerdata(:,2:end); 
cancerdata = cancerdata(ind:end,:);