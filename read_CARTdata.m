%% code to read excel files 

% INPUT: receptor = L,M,H,VH, CARTnum = 0,5,10,20 

%cancerdata = xlsread( 'Run1_HT1080.xlsx', 'Sheet1' ); 
cancerdata = xlsread( 'Run8_PBT138.xlsx', 'Sheet1' ); 


% Time in first column (in sec) 
time = cancerdata(:,1); 
time = time/3600 /24;
ind = find(time>1.75,1);
time = time(ind:end);



cancerdata  = cancerdata(:,2:end); 
cancerdata = cancerdata(ind:end,:);
%ind = find(time>80/24,1);
%time = time(1:ind);
%cancerdata = cancerdata(1:ind,:);
%time = time-time(1);
% Last cancer data is often negative. Scale that to non-negative number... 
for n = 1:size(cancerdata, 2) 
    m = size(cancerdata,1);
    while (cancerdata(m,n)<0)
        cancerdata(:,n) = cancerdata(:,n) - cancerdata(m,n); 
        m=m-1;
    end 
end



% 1 = L / 2 = M / 3 = H / 4 = VH 
receptorind = [1,1,1, 2,2,2, 3,3,3, 4,4,4 ];
receptorind = repmat(receptorind,1,4);

if( strcmp( receptor, 'L' ) )
    rind = 1; 
elseif( strcmp( receptor, 'M' ) )
    rind = 2; 
elseif( strcmp( receptor, 'H' ) )
    rind = 3; 
elseif( strcmp( receptor, 'VH' ) )
    rind = 4;
else
    disp('receptor not valid' ) 
end 

cartinit= [repmat(0,1,12),repmat(5,1,12),repmat(10,1,12),repmat(20,1,12)]; %,repmat(20,1,12),repmat(10,1,12),repmat(5,1,12),repmat(0,1,12)]; %#ok<REPMAT>


% choose data column, receptor = L,M,H,VH, CARTnum = 0,5,10,20 
cancercol = intersect( find( receptorind==rind ), find( cartinit == CARTnum ) ); 

selected_cancerdata = cancerdata(:, cancercol); 
%selected_cartdata = 



% Just take part of data 


%ind = intersect( 
%n = 2; 
%for m = 1:4 
%   cancerdata(:, (12*(m-1)+(n-1)*3+1):(12*(m-1)+n*3) ) = cancerdata(:, (12*(m-1)+(n-1)*3+1):(12*(m-1)+n*3) )*3.1616; 
%end 

