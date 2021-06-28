%clear all
%
addpath( './DRAM_Code/'); 

%% enter data to fit 
%load('PBT030-Tumor')
fitData.xdata = time - time(1); 
%C_E7_105 = [29.329;48.0114;76.6798;106.8994;136.1111;193.4253;266.1067];
%convert_size_num = (75/2 *10^(-3))^2 *pi/50;
%C_E7_105 = C_E7_105 /convert_size_num;

fitData.ydata = data;  

%% select fixed parameters 
global b %n d p g m
%a = 0.1105; 
%b = 4.2e-9; 

b = 1/6.6817;  %b=1/6 for PBT138

%% define parameters to estimate 
 params = {
    {'a', 0.4105,       0,     4} 
%     {'b', 0.3, 0.1, 10}% ParamName, starting value, uniform prior bounds     
   % {'p', 0.2,          0.1,       0.4}  
    %{'m', 3*0.1^10,   3*0.1^11,  3*0.1^9}
    %{'n', 1*0.1^7,      0.1^8,     0.1^6}
    %{'d', 0.02,         0.02,       0.06}
    %{'g', 2*0.1^7,    2*0.1^8,   2*0.1^6}
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
% figure; mcmcplot(chain,[],results,'pairs');

%%%%% plot this to check pair of chain distribution 
figure(300); mcmcplot(chain,[],results,'denspanel');

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
ind = ind(1);
params = chain(ind,:) %These are your fitted parameter values

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
[t,modFit] = ode23(@(t,y)tumor_only(t,y,params), plottime, fitData.ydata(1,:)');


figure;
plot(plottime,modFit,'-r','LineWidth',1)
hold on
plot(fitData.xdata,fitData.ydata(:,1),'ok','MarkerSize',6,'MarkerFaceColor','k')

xlabel('Time','FontSize',14)
ylabel('Tumor Size','FontSize',14)
set(gca,'FontSize',14)


function SS = LLHfunc(params,fitData)

[t,tumVol] = ode23(@(t,y)tumor_only(t,y,params), fitData.xdata, fitData.ydata(1,:)' );

%% Choose 'Sum of squares' 
% SS = sum((tumVol(:,1) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
nsig = 0.1;  %%%% nsig*100 percent noise 
sigma =  fitData.ydata(:,1) * nsig; 
SS = -sum(log((1./(sqrt(2*pi)*sigma)).*exp(-((tumVol(:,1) - fitData.ydata(:,1)).^2)./(2*sigma.^2))));


end