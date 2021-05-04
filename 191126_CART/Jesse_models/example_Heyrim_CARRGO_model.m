
%% Test MCMC

%initParams
%fit_MCMC_DRAM(obj, data, initParams, options)

%%%% The sum-of-squares "ssfun" method in the model class needs "y0" as a
%%%% global variable.
%%%% Note that I made RHS_CARRGO so that CAR T cell dynamics are "paused",
%%%% so I can just include them in initial conditions
global y0; y0 = [1 1];

%%%% The onset of CAR T cell dyamics (including killing) is controlled for
%%%% t < t0_CART
global t0_CART; t0_CART = 1;


%%%% simulated noise levels (CANCER, CAR T) for observations
noise_sigma = [0.01 0.05];

%%%% create the model object
modelName = "CARRGO";
myModel = model(modelName);

%%%% the "true" parameter values for kappa_1, kappa_2, theta
%%%% we use these to create the simulated data
pvect = [3.5 2 0.02];

%%%% create simulated data (true and noisy observations)
data.xdata = 0:0.1:5;
data.x{1} = true([1 length(data.xdata)]);
data.x{2} = false([1 length(data.xdata)]);
data.ydata = zeros([length(data.xdata) 2]);

dd = data;
dd.x{1}(:) = 1;
dd.x{2}(:) = 1;

y_true_observed = myModel.modelfun(data, pvect);
y_true_full = myModel.modelfun(dd, pvect);

noise = randn(size(y_true_observed))*diag(noise_sigma);
% y_noisy = y_true + noise;
y_noisy_observed = y_true_observed.*(1+noise);
y_noisy_full = y_true_full.*(1+noise);

data.ydata = y_noisy_observed;


%%%% pick options

startpvect = [4 3.3 .04];
[results, chain, s2chain, sschain] = myModel.fit_MCMC_DRAM(data, startpvect);

theta = results.theta;
y_fit = myModel.modelfun(dd, theta);


figure(4)
hold on;
x = data.xdata;
for(jj = 1:2)
    subplot(2,1,jj)
    plot(x,y_true_full(:,jj),x,y_noisy_full(:,jj), x,y_fit(:,jj))
    legend(["True Orbit", "Noisy Orbit", "Fit Orbit"])
end


hold off