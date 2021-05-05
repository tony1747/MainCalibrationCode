%% code to read excel files 

%cancerdata = xlsread( 'Run1_HT1080.xlsx', 'Sheet1' ); 
cancerdata = xlsread( 'Run8_PBT138.xlsx', 'Sheet1' ); 
receptor = {'L','L','L','M','M','M','H','H','H','VH','VH','VH'};
receptor = repmat(receptor,1,8);
cartinit= [repmat(0,1,12),repmat(5,1,12),repmat(10,1,12),repmat(20,1,12),repmat(20,1,12),repmat(10,1,12),repmat(5,1,12),repmat(0,1,12)]; %#ok<REPMAT>

% Time in first column (in sec) 
time = cancerdata(:,1); 
time = time/3600 /24;
ind = find(time>1.25,1);
time = time(ind:end);

cancerdata  = cancerdata(:,2:end); 
cancerdata = cancerdata(ind:end,:);