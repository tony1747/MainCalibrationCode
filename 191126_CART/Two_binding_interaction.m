%clear all
%
addpath( './DRAM_Code/'); 

%% enter data to fit 
% read_CARTdata; 
CancerType = 'L'; 
global CARTratio
CARTratio = 0.2; 

[time, data] = load_CARTdata( CancerType, CARTratio ); 

fitData.xdata = time - time(1);
fitData.ydata(:,1) = data(:,1);  

%% select fixed parameters 
global a b %n d p g m
[a,b] = set_CancerGrowthParams( CancerType ); 

%% define parameters to estimate 
 params = {
   % {'a', 0.4105,       0.168,     4} 
  %  {'b', 0.3, 0.1, 10}% ParamName, starting value, uniform prior bounds     
    {'p', 0.4862,      0,       5}  
    %{'m', 0.1,     0,       0.5}
    %{'n', 1 ,     0,        5}
    {'d', 0.0598,     0,       0.1}
    {'g', 14.0247,      1,   20}
    {'alpha', 2, 1, 20}
    {'beta', 0.025, 0, 1}
    {'gamma', 0.001, 0, 0.5}
    {'delta', 1, 0, 5}
    {'epsilon', 27.85, 0, 60}
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
%if( length( unique(chain(:)) ) > 8 ) 
    % figure(301); mcmcplot(chain,[],results,'denspanel');
%end

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
ind = ind(1);
params = chain(ind,:) %These are your fitted parameter values

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
%plottime = [0:.1:30]; 

[t,modFit] = ode23(@(t,y)two_binding(t,y,params), plottime, [fitData.ydata(1,1)*CARTratio,fitData.ydata(1,1)]');


figure; hold on; 
plot(plottime,modFit(:,1),'--b','LineWidth',1)
plot(plottime,modFit(:,2),'-r','LineWidth',1)
plot(fitData.xdata,fitData.ydata(:,1),'ok','MarkerSize',6,'MarkerFaceColor','k')
error=norm(modFit-fitData.ydata(:,1));
xlabel('Time','FontSize',14)
ylabel('Tumor Size','FontSize',14)
set(gca,'FontSize',14)


function SS = LLHfunc(params,fitData)
global CARTratio
[t,tumVol] = ode23(@(t,y)two_binding(t,y,params), fitData.xdata, [fitData.ydata(1,1)*CARTratio,fitData.ydata(1,1)]' );

%% Choose 'Sum of squares' 
% SS = sum((tumVol(:,1) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
%nsig = 0.1;  %%%% nsig*100 percent noise 
%sigma =  fitData.ydata(:,1) * nsig; 
SS = norm((tumVol(:,2) - fitData.ydata(:,1)));


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
