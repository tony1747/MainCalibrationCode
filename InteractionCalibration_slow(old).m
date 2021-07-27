%clear all
%
addpath( './DRAM_Code/'); 
%% select fixed parameters 
global a b CARTratio d %n d p g m
%% enter data to fit 
%load('PBT030-Tumor')
fitData.xdata = time;
%C_E7_105 = [29.329;48.0114;76.6798;106.8994;136.1111;193.4253;266.1067];
%convert_size_num = (75/2 *10^(-3))^2 *pi/50;
%C_E7_105 = C_E7_105 /convert_size_num;

fitData.ydata(:,1)= data;
fitData.ydata(1,2)= data(1)*CARTratio;
fitData.ydata(1,3)=0;

[a,b] = set_CancerGrowthParams( receptor );
d=0.0412;



%% define parameters to estimate 
 params = {
   % {'a', 0.4105,       0.168,     4} 
    
    %{'b', 0.3, 0.1, 10}% ParamName, starting value, uniform prior bounds     
    %{'p', 0.2,          0.1,       0.4}  
    {'k_1',0.3,0,10}
    {'alpha',0.3,0,20}
    {'beta',0.3,0,10} %alpha=k_n1+k_3
    {'gamma',0.3,0,10}  %beta=k_n1+k_2
    %{'k_n1',0.1,0,1}  %gamma=k_n1+k_2+k_3
    {'p', 0.4, 0, 3}
    %{'m', 3*0.1^10,   3*0.1^11,  3*0.1^9}
    %{'m', 0, 0, 0.1^14}
    %{'n', 1*0.1^7,      0.1^8,     0.1^6}
    %{'n', 0, 0, 0.1^14}
    %{'d', 0.02,         0.02,       0.06}
    %{'d', 0, 0, 0.1^14}
    %{'g', 2*0.1^7,    2*0.1^8,   2*0.1^6}
    {'g', 2, 1, 20}
    };

%% This is the likelihood function 
model.ssfun = @LLHfunc; 
% SSQfunc
options.nsimu = 10000; %number of samples for DRAM

%% Run mcmc 
[results, chain, s2chain, ss2chain] = mcmcrun(model,fitData,params,options);



%% Plot results 

%%%%% plot to check chain convergence 
% figure; mcmcplot(chain,[],results,'chainpanel');

%%%%% plot to check pair of chain relation  
figure; mcmcplot(chain,[],results,'pairs');

%%%%% plot this to check pair of chain distribution 
%figure; mcmcplot(chain,[],results,'denspanel');

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
error = ss2chain(ind(1));
ind = ind(1);
params = chain(ind,:); %These are your fitted parameter values
slow_1_binding=strcat('Slow_1_binding',receptor,int2str(CARTnum),'_',int2str(n),'.mat');

save(slow_1_binding,'params','error')

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
%[t,modFit] = ode23(@(t,y)tumor_cart(t,y,params), plottime, [fitData.ydata(1,1)*0.2,fitData.ydata(1,1)]');
[t,modFit] = ode23(@(t,y)tumor_cart(t,y,params), plottime, fitData.ydata(1,:)');


figure;hold on;
plot(plottime,modFit(:,1)+modFit(:,3),'--r','LineWidth',1)
plot(plottime,modFit(:,2)+modFit(:,3),'-b','LineWidth',1)
%plot(plottime,modFit(:,3),'
plot(fitData.xdata,fitData.ydata(:,1),'ok','MarkerSize',6,'MarkerFaceColor','k')
legend({'Cancer','CAR-T','data'})

xlabel('Time','FontSize',14)
ylabel('Tumor Size','FontSize',14)
set(gca,'FontSize',14)
slow_1_binding=strcat('Slow_1_binding',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
saveas(gcf,slow_1_binding)
close all

function SS = LLHfunc(params,fitData)

%[t,tumVol] = ode23(@(t,y)tumor_cart(t,y,params), fitData.xdata, [fitData.ydata(1,1)*0.2,fitData.ydata(1,1)]' );
[t,tumVol] = ode23(@(t,y)tumor_cart(t,y,params), fitData.xdata, fitData.ydata(1,:)' );

%% Choose 'Sum of squares' 
SS = sum((tumVol(:,1)+tumVol(:,3) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
%nsig = 0.1;  %%%% nsig*100 percent noise 
%sigma =  fitData.ydata(:,1) * nsig; 
%SS = -sum(log((1./(sqrt(2*pi)*sigma)).*exp(-((tumVol(:,1) - fitData.ydata(:,1)).^2)./(2*sigma.^2))));


end

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

function [a,b]=set_CancerGrowthParamsHT( CancerType )

switch (CancerType )
    case 'L'
        a=2.6136;
        b=0.1660;
    case 'M'
        a=2.1991;
        b=0.2861;
    case 'H'
        a=2.3474;
        b=0.3247;
    case 'VH'
        a=1.8704;
        b=0.4702;
end
end
