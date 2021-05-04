%% code to read excel files 

% cancerdata = xlsread( 'Run1_HT1080.xlsx', 'Sheet1' ); 
cancerdata = xlsread( 'Run8_PBT138.xlsx', 'Sheet1' ); 


% Time in first column (in sec) 
time = cancerdata(:,1); 
time = time/3600 /24; 

cancerdata  = cancerdata(:,2:end); 

% Last cancer data is often negative. Scale that to non-negative number... 
for n = 1:size(cancerdata, 2) 
    if(cancerdata(end,n)<0)
        cancerdata(:,n) = cancerdata(:,n) - cancerdata(end,n); 
    end 
end

% Select time > 1.75 day part of data 
tind = find( time > 1.25, 1 ); 
time = time(tind:end); 
cancerdata = cancerdata(tind:end, :); 



% Just take part of data 


%ind = intersect( 
%n = 2; 
%for m = 1:4 
%   cancerdata(:, (12*(m-1)+(n-1)*3+1):(12*(m-1)+n*3) ) = cancerdata(:, (12*(m-1)+(n-1)*3+1):(12*(m-1)+n*3) )*3.1616; 
%end 

