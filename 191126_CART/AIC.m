params=xlsread( 'fitted params (PBT138, slow, k_n1 included).xlsx', 'Sheet1' );
global a b d
AIC_list=[];
for receptor = ['L', 'M', 'H']
    for CARTnum = [5, 10, 20]
        read_CARTdata() %% <- change  last part, data to cancerdata <- maybe you can also make something to select the cancer cell type outside.
        CARTratio=1/CARTnum;
        for n = 1:3
            data = selected_cancerdata(:,n);
            if receptor=='M'
                data=data*3;
            end
            fitData.xdata = time;
            %C_E7_105 = [29.329;48.0114;76.6798;106.8994;136.1111;193.4253;266.1067];
            %convert_size_num = (75/2 *10^(-3))^2 *pi/50;
            %C_E7_105 = C_E7_105 /convert_size_num;
            
            fitData.ydata(:,1)= data;
            fitData.ydata(1,2)= data(1)*CARTratio;
            %fitData.ydata(1,3)=0;
            
            [a,b] = set_CancerGrowthParams( receptor );
            d=0.0412;
            if receptor=='L' & CARTnum==5
                paramsdata = params(1:3,:);
            elseif receptor=='L' & CARTnum==10
                paramsdata = params(4:6,:);
            elseif receptor=='L' & CARTnum==20
                paramsdata = params(7:9,:);
            elseif receptor=='M' & CARTnum==5
                paramsdata = params(10:12,:);
            elseif receptor=='M' & CARTnum==10
                paramsdata = params(13:15,:);
            elseif receptor=='M' & CARTnum==20
                paramsdata = params(16:18,:);
            elseif receptor=='H' & CARTnum==5
                paramsdata = params(19:21,:);
            elseif receptor=='H' & CARTnum==10
                paramsdata = params(22:24,:);
            else 
                paramsdata = params(25:27,:);
            end
            paramsdata_n = params(n,:);
            [t,tumVol] = ode23(@(t,y)tumor_cart_only(t,y,paramsdata_n), fitData.xdata, fitData.ydata(1,:)' );
            chi_sq = sum(((fitData.ydata(:,1)- tumVol(:,1))./tumVol(:,1)).^2);
            Akaike = length(tumVol)*log(chi_sq/length(tumVol))+2*(7+1); %use 7+1 for fast reaction
            %use 9+1 for slow reaction
            AIC_list=[AIC_list; Akaike];
        end
        
        %% here save parameters/posterior
    end
end
csvwrite("AIC_list (fast_binding).csv",AIC_list);
function [a,b] = set_CancerGrowthParams( CancerType )

switch( CancerType ) 
    case 'L'
        a = 0.3953;          
    case 'M'
        a = 0.3606;  
    case 'H'
        a = 0.2187;  
end

b = 1/6; 

end 