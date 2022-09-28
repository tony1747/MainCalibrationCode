%clear all
%
addpath( './DRAM_Code/'); 

%% enter data to fit 
load( 'Run8_PBT138_data.mat' ); 
fitData.xdata = time;
data = cancerdata(:,[13,25,37]); 

fitData.ydata(:,1) = data(:,1);
fitData.ydata(:,2) = data(:,2);
fitData.ydata(:,3) = data(:,3);

%% select fixed parameters 
global a b %n d p g m
% Run 8 L 
a = 0.3953;  
b = 1/6; 


%% define parameters to estimate 
 params = {
   % {'a', 0.4105,       0.168,     4} 
  %  {'b', 0.3, 0.1, 10}% ParamName, starting value, uniform prior bounds     
    {'p', 0.2,      0,       1.0}  
    {'m', 0.01,     0,       0.5}
    {'n', 3 ,     0,        5}
%     {'d', 0.02,     0,       0.06}
    {'g', 10,      1,   20}
    };
global d p 
d = 0.05; 
% p = 0.2; 

%% This is the likelihood function 
model.ssfun = @LLHfunc; 
% SSQfunc
options.nsimu = 1*10^4; %number of samples for DRAM

%% Run mcmc 
[results, chain, s2chain, ss2chain] = mcmcrun(model,fitData,params,options);



%% Plot results 

%%%%% plot to check chain convergence 
figure; mcmcplot(chain,[],results,'chainpanel');

%%%%% plot to check pair of chain relation  
figure; mcmcplot(chain,[],results,'pairs');

%%%%% plot this to check pair of chain distribution 
figure; mcmcplot(chain,[],results,'denspanel');

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
ind = ind(1);
params = chain(ind,:) %These are your fitted parameter values

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
CARTratio = [0.2, 0.1, 0.05]; 

figure; 
for i = 1:3 
    [t,modFit] = ode23(@(t,y)tumor_cart_only(t,y,params), plottime, [fitData.ydata(1,i)*CARTratio(i),fitData.ydata(1,i)]');
    
    subplot(1, 3, i); hold on; 
    plot(plottime,modFit(:,1),'--b','LineWidth',1)
    plot(plottime,modFit(:,2),'-r','LineWidth',1)
    plot(fitData.xdata,fitData.ydata(:,i),'ok','MarkerSize',6,'MarkerFaceColor','k')

    xlabel('Time','FontSize',14)
    ylabel('Tumor Size','FontSize',14)
    set(gca,'FontSize',14)   
end 



function SS = LLHfunc(params,fitData)
tumVol = zeros(length(fitData.xdata),3);
[t,tumVoltmp] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,1)*0.2,fitData.ydata(1,1)]' );
tumVol(:,1) = tumVoltmp(:,2); 
[t,tumVoltmp] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,2)*0.1,fitData.ydata(1,2)]' );
tumVol(:,2) = tumVoltmp(:,2); 
[t,tumVoltmp] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,3)*0.05,fitData.ydata(1,3)]' );
tumVol(:,3) = tumVoltmp(:,2); 

%% Choose 'Sum of squares' 
% SS = sum((tumVol(:,1) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
nsig = 0.1;  %%%% nsig*100 percent noise 
sigma =  fitData.ydata * nsig; 
SS = -sum(sum(log((1./(sqrt(2*pi)*sigma)).*exp(-((tumVol - fitData.ydata).^2)./(2*sigma.^2)))));


end