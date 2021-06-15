%clear all
%
addpath( './DRAM_Code/'); 

%% enter data to fit
% read_CARTdata; 
%CancerType = 'L'; 
global CARTratio
%CARTratio = 0.2; 
%use a for loop
%for CancerType=['L','M','H']
%for CARTratio=[0.2,0.1,0.05]
%disp({CancerType,CARTratio})
%InteractionCalibration()
%end
%end
%[time, data] = load_CARTdata( CancerType, CARTratio ); 

fitData.xdata = time - time(1);
fitData.ydata(:,1) = data; 
fitData.ydata(1,2) = data(1)*CARTratio;

%% select fixed parameters 
global a b d %n p g m
[a,b] = set_CancerGrowthParams( receptor ); 
d=0.0412;

%% define parameters to estimate 
 params = {
   % {'a', 0.4105,       0.168,     4} 
  %  {'b', 0.3, 0.1, 10}% ParamName, starting value, uniform prior bounds     
    {'p', 0.4,      0,       3}  
    {'m', 0.1,     0,       0.5}
    {'n', 1 ,     0,        10}
    %{'d', 0.02,     0,       0.5}
    {'g', 2,      1,   20}
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
figure(300); mcmcplot(chain,[],results,'pairs');

%%%%% plot this to check pair of chain distribution 
if( length( unique(chain(:)) ) > 4 )  %this number is dependent on the number of parameters we want to calibrate
    figure(301); mcmcplot(chain,[],results,'denspanel');
end

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
error = ss2chain(ind(1));
ind = ind(1);
params = chain(ind,:); %These are your fitted parameter values
fast_1_binding=strcat('Fast_1_binding',receptor,int2str(CARTnum),'_',int2str(n),'.mat');
save(fast_1_binding,'params','error')
if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 )
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
%plottime = [0:.1:30]; 

[t,modFit] = ode23(@(t,y)tumor_cart_only(t,y,params), plottime, fitData.ydata(1,:)');


figure;  hold on;
plot(plottime,modFit(:,1),'--r','LineWidth',1)
plot(plottime,modFit(:,2),'-b','LineWidth',1)
plot(fitData.xdata,fitData.ydata(:,1),'ok','MarkerSize',6,'MarkerFaceColor','k')
legend({'Cancer','CAR-T','data'})

xlabel('Time','FontSize',14)
ylabel('Tumor Size','FontSize',14)
set(gca,'FontSize',14)
fast_1_binding=strcat('Fast_1_binding',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
saveas(gcf,fast_1_binding)
close all
%fi
function SS = LLHfunc(params,fitData)
global CARTratio
[t,tumVol] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, fitData.ydata(1,:)' );

%% Choose 'Sum of squares' 
% SS = sum((tumVol(:,1) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
%nsig = 0.1;  %%%% nsig*100 percent noise 
%sigma =  fitData.ydata(:,1) * nsig; 
%SS = -sum(log((1./(sqrt(2*pi)*sigma)).*exp(-((tumVol(:,2) - fitData.ydata(:,1)).^2)./(2*sigma.^2))));
SS = norm((tumVol(:,1) - fitData.ydata(:,1)));

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


