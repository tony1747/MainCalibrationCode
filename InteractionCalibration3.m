%clear all
%
addpath( './DRAM_Code/'); 

%% enter data to fit 
%load('PBT030-Tumor')
fitData.xdata = time;
%C_E7_105 = [29.329;48.0114;76.6798;106.8994;136.1111;193.4253;266.1067];
%convert_size_num = (75/2 *10^(-3))^2 *pi/50;
%C_E7_105 = C_E7_105 /convert_size_num;

fitData.ydata(:,1) = data(:,1);
fitData.ydata(:,2) = data(:,2);
fitData.ydata(:,3) = data(:,3);

%% select fixed parameters 
global a b %n d p g m
a = 2.5874; 
b = 0.4001; 

%% define parameters to estimate 
  params = {
   % {'a', 0.02,       0.01,     0.52} 
    
    % {'b', 1e-7, 1e-14, 1e-4}% ParamName, starting value, uniform prior bounds     
    %{'p', 0.2,          0.1,       0.4}  
    {'p', 0.2, 0.1, 0.4}
    %{'m', 3*0.1^10,   3*0.1^11,  3*0.1^9}
    {'m', 0.2, 0, 0.5}
    %{'n', 1*0.1^7,      0.1^8,     0.1^6}
    {'n', 2, 0, 10}
    %{'d', 0.02,         0.02,       0.06}
    {'d', 0.005, 0, 0.01}
    %{'g', 2*0.1^7,    2*0.1^8,   2*0.1^6}
    {'g', 1.5, 1, 20}
    };

%% This is the likelihood function 
model.ssfun = @LLHfunc; 
% SSQfunc
options.nsimu = 2*10^4; %number of samples for DRAM

%% Run mcmc 
[results, chain, s2chain, ss2chain] = mcmcrun(model,fitData,params,options);



%% Plot results 

%%%%% plot to check chain convergence 
% figure; mcmcplot(chain,[],results,'chainpanel');

%%%%% plot to check pair of chain relation  
figure; mcmcplot(chain,[],results,'pairs');

%%%%% plot this to check pair of chain distribution 
figure; mcmcplot(chain,[],results,'denspanel');

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain));
ind = ind(1);
params = chain(ind,:) %These are your fitted parameter values
noCAR_T=strcat('noCAR-T',int2str(n_cell),'.mat');
save(noCAR_T,'params')

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 

plottime = fitData.xdata; 
[t,modFit] = ode23(@(t,y)tumor_cart_only(t,y,params), plottime, [fitData.ydata(1,1)*0.2,fitData.ydata(1,1)]');


figure;
plot(plottime,modFit,'-r','LineWidth',1)
hold on
plot(fitData.xdata,fitData.ydata(:,1),'ok','MarkerSize',6,'MarkerFaceColor','k')

xlabel('Time','FontSize',14)
ylabel('Tumor Size','FontSize',14)
set(gca,'FontSize',14)


function SS = LLHfunc(params,fitData)
tumVol = zeros(length(fitData.xdata),3);
[t,tumVol(:,1)] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,1)*0.2,fitData.ydata(1,1)]' );
[t,tumVol(:,2)] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,2)*0.1,fitData.ydata(1,2)]' );
[t,tumVol(:,3)] = ode23(@(t,y)tumor_cart_only(t,y,params), fitData.xdata, [fitData.ydata(1,3)*0.05,fitData.ydata(1,3)]' );

%% Choose 'Sum of squares' 
% SS = sum((tumVol(:,1) - fitData.ydata(:,1)).^2);

%% Choose Gaussian noise of sigma variance 
nsig = 0.1;  %%%% nsig*100 percent noise 
sigma =  fitData.ydata * nsig; 
SS = -sum(sum(log((1./(sqrt(2*pi)*sigma)).*exp(-((tumVol - fitData.ydata).^2)./(2*sigma.^2)))));


end