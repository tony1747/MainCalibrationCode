function [time, data] = load_CARTdata( receptorDensity, nCART0 ) 

load( 'Run8_PBT138_data(longer time)' ); 

switch(receptorDensity) 
    case 'L'
        ctype = 1; 
    case 'M'
        ctype = 2; 
    case 'H'
        ctype = 3; 
end

ind = intersect( find(dataCANCERnum == ctype), find(dataCARTnum == nCART0) ); 
ind = ind( ind < 49 ); %take front half 

data = cancerdata(:,ind); 
data = data(:,1);
data = data+data(end-1);
end 

