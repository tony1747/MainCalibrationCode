


clear all;

%% Test: model implement with "CARRGO" model, and generate, plot output
global y0; y0 = [1 1];
global t0_CART; t0_CART = 1;

data.xdata = 0:0.1:5;
data.x{1} = true([1 length(data.xdata)]);
data.x{2} = true([1 length(data.xdata)]);
data.ydata = zeros([length(data.xdata) 2]);


modelName = "CARRGO";

myModel = model(modelName);

if(~exist("pvect",'var'))
    pvect = myModel.par_default(myModel.par_inds);
end

outy = myModel.modelfun(data, pvect);


%%% Plot it
figure(1);

for ii=1:2
    subplot(1,2,ii)
    plot(data.xdata,outy(:,ii))
end



%% Now we test misc functions associated with MCMC

pvect = [5.5 5.5];
inds = [3 5];

paramsStruct1 = myModel.get_paramsForMCMC();

disp(" We expect an error here, due to improper usage.")
try
    paramsStruct2 = myModel.get_paramsForMCMC(pvect);
catch me
    disp( getReport( me, 'basic', 'hyperlinks', 'on' ) )
end

paramsStruct3 = myModel.get_paramsForMCMC(pvect,inds);

pvect = [5.5 5.5 0.0004];

disp(" We fix the error here, by matching pvect to provided default indices.")
paramsStruct2 = myModel.get_paramsForMCMC(pvect);



%%% Test Sum of Squares function


global y0; y0 = [1 1];
global t0_CART; t0_CART = 1;

noise_Sigma = 0:0.01:0.2;
noise_samples = 3;
plotPairs_ss_sigma = zeros([length(noise_Sigma),1+noise_samples]);

modelName = "CARRGO";
myModel = model(modelName);
if(~exist("pvect",'var'))
    pvect = myModel.par_default(myModel.par_inds);
end

data.xdata = 0:0.1:5;
data.x{1} = true([1 length(data.xdata)]);
data.x{2} = true([1 length(data.xdata)]);
data.ydata = zeros([length(data.xdata) 2]);

y_true = myModel.modelfun(data, pvect);


for ii=1:length(noise_Sigma)
   sigma = noise_Sigma(ii);
   for jj = 1:noise_samples
       noise = randn(size(y_true))*sigma;
       y_noisy = y_true + noise;
       data.ydata = y_noisy;
       ss = myModel.ssfun(data, pvect, [1]);
       
       plotPairs_ss_sigma(ii,1) = sigma;
       plotPairs_ss_sigma(ii,1+jj) = ss;
   end
end

figure(2);
clf;
hold on;
for jj = 1:noise_samples
    plot(plotPairs_ss_sigma(:,1),plotPairs_ss_sigma(:,1+jj))
end
legend(strcat("Sample ", string(1:noise_samples)))
title("SS Function for noisy orbit Y_hat")
xlabel("Noise Level (sigma)")
ylabel("SSfun")
hold off


%% fit fmincon



global y0; y0 = [1 1];
global t0_CART; t0_CART = 1;

%noise_Sigma = 0:0.01:0.2;
%noise_samples = 3;
noise_sigma = [0.01 0.05];
%plotPairs_ss_sigma = zeros([length(noise_Sigma),1+noise_samples]);

modelName = "CARRGO";
myModel = model(modelName);
if(~exist("pvect",'var'))
    pvect = myModel.par_default(myModel.par_inds);
end

data.xdata = 0:0.1:5;
data.x{1} = true([1 length(data.xdata)]);
data.x{2} = true([1 length(data.xdata)]);
data.ydata = zeros([length(data.xdata) 2]);

y_true = myModel.modelfun(data, pvect);

noise = randn(size(y_true))*diag(sigma);
% y_noisy = y_true + noise;
y_noisy = y_true.*(1+noise);

data.ydata = y_noisy;

[myfit ssout] = myModel.fit_fmincon(data, 30);


%% Test MCMC

%initParams
%fit_MCMC_DRAM(obj, data, initParams, options)


global y0; y0 = [1 1];
global t0_CART; t0_CART = 1;


noise_sigma = [0.01 0.05];


modelName = "CARRGO";
myModel = model(modelName);

pvect = [3.5 2 0.02];


data.xdata = 0:0.1:5;
data.x{1} = true([1 length(data.xdata)]);
data.x{2} = true([1 length(data.xdata)]);
data.ydata = zeros([length(data.xdata) 2]);

y_true = myModel.modelfun(data, pvect);

noise = randn(size(y_true))*diag(noise_sigma);
% y_noisy = y_true + noise;
y_noisy = y_true.*(1+noise);

data.ydata = y_noisy;

%options.nsimu = 1e4;
startpvect = [4 3.3 .04];
[results, chain, s2chain, sschain] = myModel.fit_MCMC_DRAM(data, startpvect);

theta = results.theta;
y_fit = myModel.modelfun(data, theta);


figure(4)
hold on;
x = data.xdata;
for(jj = 1:2)
    subplot(2,1,jj)
    plot(x,y_true(:,jj),x,y_noisy(:,jj), x,y_fit(:,jj))
    legend(["True Orbit", "Noisy Orbit", "Fit Orbit"])
end


hold off
