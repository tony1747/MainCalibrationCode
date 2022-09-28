classdef model
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model_name string %% every model should have a string NAME, associated with files: paramInfo_NAME.csv and RHS_NAME.m
        par_num {mustBeInteger} = 5; %% number of parameters
        par_names string %% names of parameters (ORDER MATTERS)
        par_ub double %% allowed upper bounds of parameters
        par_lb double  %% allowed lower bounds of parameters
        par_default double %% default values of parameters. Crucial for fixed parameters
        par_inds {mustBeInteger}; %% default indices of non-fixed parameters
        par_struct cell %% this is a the parameter struct object needed for MCMC
        y_num {mustBeInteger} = 2; %% number of state variables
        y_obs {mustBeInteger} = [1]; %% indices of observed state variables
        RHS_handle function_handle %% function for RHS of ODE
    end
    
    methods
        
        %%% constructor, that calls other methods
        function obj = model(modelName)
            if(nargin == 1)
            obj = model();
            obj = obj.model_fromFile(modelName);
            end
            %%% if no argument, makes a blank model
        end
        
        
        %%%%%% finished, works (4/13/20)
        function obj = model_fromFile(obj, modelName)

            obj.model_name = modelName;
            
            RHSfile_prefix = strcat("RHS_",modelName);
            parfn_prefix = strcat("paramInfo_",modelName);
            
            obj.RHS_handle = str2func(RHSfile_prefix);
            
%             testfn =str2func(RHSfile_prefix);
            

%             %%%% This bit here was when I saved parameter info
%             %%%% as a function, instead of csv. Obsolete
%             parfn = str2func(parfn_prefix);
%             parInfo = parfn();
            parInfo = paramInfo(modelName);
            
            

            N_p = length(parInfo.names);
            
            
            %%% get the parameter info
            obj.par_names = parInfo.names;
            obj.par_num = N_p;
            obj.par_ub = parInfo.ub;
            obj.par_lb = parInfo.lb;
            
            if(isfield(parInfo,"default"))
                obj.par_default = parInfo.default;
            else
                obj.par_default = (parInfo.ub+parInfo.lb)/2;
            end
            
            if(isfield(parInfo,"inds"))
                    obj.par_inds = parInfo.inds;
                else if (isfield(parInfo,"pind"))
                    obj.par_inds = parInfo.pind;
                else
                    obj.par_inds = (3:N_p);
                end
            end
            
            if(isfield(parInfo,"y_num"))
                
            end
            
            %%% make structured cell of parameters
            obj = obj.make_parStruct();

            
        end
        
        %%%%% parameter related functions
        
        %%%%%% done, untested
        function obj = make_parStruct(obj)
            %%%% This function makes the par_struct property
            %%%% from other properties.  par_struct is used
            %%%% below in "get_paramsForMCMC"

            
            names = obj.par_names;
            ub = obj.par_ub;
            lb = obj.par_lb;
            sv = obj.par_default;
            N_p = length(names);

            tmp = cell(repmat({cell(1,4)},[N_p,1]));

            for ii=1:length(tmp)
                tmp{ii} = {names(ii), sv(ii), lb(ii),ub(ii)};
            end
            
            obj.par_struct = tmp;
        end
        
        
        %%%%%% done, untested?
        %%%%%% ERROR, see test_modelclass file, figure out
        function params = get_paramsForMCMC(obj, pvect, ind)
            %%%% This function makes a cell structure
            %%%% for parameters and their bounds,
            %%%% for use in MCMC
            
            bool_ind = exist('ind','var');
             if ~bool_ind
               ind = obj.par_inds;
             end
             
             if ~exist('pvect','var')
               pvect = obj.par_default(ind);
             else
                 if (length(pvect) ~= length(ind))
                     msg = strcat( string(length(pvect)), " parameters provided without indices.", " Default # of indices: ", string(length(ind)) );
                     error(msg);
                 end
             end
            
             
             
             tmp = obj.par_struct;
             

             
             params = cell(length(ind),1); 
             for n = 1:length(ind) 
                tmp{ind(n)}{2} = pvect(n); 
                params{n} = tmp{ind(n)};  
             end
                
        end
        
        function parvect = get_params(obj, pvect, ind)
             if ~exist('ind','var')
               ind = obj.par_inds;
             end
             
            if (length(pvect) ~= length(ind))
                msg = strcat("Error. Method parvect given ",string(length(pvect)),...
                    " parameters, expected ",string(length(ind)));
               error(msg);
            end
            
             parvect = obj.par_default;
             parvect(ind) = pvect;
            
        end
        
        
        %%% Useful functions: modelfun, ssfun
       
        function varargout = modelfun(obj, data, params, resolve)
                % this function creates trajectories,
                % given time info from "data" struct
                % and (non-fixed) parameters "params"
                % Note that pind/pfix depends on "obj" attribute
                global y0
                
%                 disp("Params")
%                 disp(params)
%                 disp("DATA")
%                 disp(data)

                
                M = size(data.ydata,1);
                %%% insert a check here
                
                 if ~exist('params','var')
                   tmpind = obj.par_inds;
                   params = obj.par_default(tmpind);
                 end
                 
                 if ~exist('resolve','var')
                     resolve = "none";
                 end

                 
                
                parvect = obj.get_params(params);
                tspan = data.xdata;
            
                sol = ode45(@(t,y) obj.RHS_handle(t,y, parvect), tspan,  transpose(y0));
                ysol = zeros( size( data.ydata ) );
                
                no_resolve = false;
                switch resolve
                    case {'union', 'or'}
                        disp(["case 1" ; resolve])
                        tmp = reshape([data.x{:}],size(data.ydata));
                        inds = arrayfun(@(ii) any(tmp(ii,:)), 1:M);
                    case {'intersect','intersection'}
                        disp(["case 2" ; resolve])
                        tmp = reshape([data.x{:}],size(data.ydata));
                        inds = arrayfun(@(ii) all(tmp(ii,:)), 1:M);
                    case {'max'}
                        disp(["case 3" ; resolve])
                        inds = true(size(data.x{1}));
                    otherwise
                        inds = "none";
                        no_resolve = true;
                end
                
                if(~no_resolve)
                    for ii=1:obj.y_num
                        ysol(inds,ii) = interp1( sol.x, sol.y(ii,:), data.xdata(inds) );
                    end
                else
                    for ii=1:obj.y_num
                        ysol(data.x{ii},ii) = interp1( sol.x, sol.y(ii,:), data.xdata(data.x{ii}) );
                    end
                end
                %ysol(data.x{1},2) = interp1( sol.x, sol.y(2,:), data.xdata(data.x{1}) );
                
                %y_num
            
                varargout{1} = ysol;
                varargout{2} = sol;
        end
        
        
        %%% make this into an SS function
        function varargout = ssfun(obj, data, params, observed_variables)
                % this function finds,
                % given time info from "data" struct
                % and (non-fixed) parameters "params"
                % Note that pind/pfix depends on "obj" attribute
                global y0

                
                 if ~exist('observed_variables','var')
                   observed_variables = obj.y_obs;
                 end

                ysol = obj.modelfun(data, params, "union");
                
                ydata = data.ydata(:,observed_variables);
                
                boolmatrix = reshape([data.x{:}], size(data.ydata));
                indmatrix = false(size(boolmatrix));
                
                indmatrix(:,observed_variables) = boolmatrix(:,observed_variables);
           
                %%% code from ssfun
                %%% function ss = SSQfunc_CARRGO(params, data)
                %soly = modelfun_CARRGO_modified(data, params);
                %ss = sum(temp(data.x{1}).^2);
                %%% end
                
                temp = ysol(indmatrix) - data.ydata(indmatrix);
                ss = sum(temp.^2);
                
                
                varargout{1} = ss;
                varargout{2} = temp;
                
                
        end
        
        

            
            
        function out = get_LHS(obj, paraminput)
            
            if(length(paraminput) == 1 & isnumeric(paraminput))
                N_samples = paraminput;
            end
            
            N_par_vary = length(obj.par_inds);
            disp([N_samples N_par_vary])
            lhs = lhsdesign(N_samples, N_par_vary);
            
            lb = obj.par_lb(obj.par_inds);
            ub = obj.par_ub(obj.par_inds);
            range = ub-lb;
            lhs_params = repmat(lb, [N_samples 1]) + lhs*diag(range);
            out = lhs_params;
            
        end
            

        function varargout = fit_fmincon(obj, data, paraminput)
            
            if(length(paraminput) == 1 & isnumeric(paraminput))
                LHS = obj.get_LHS(paraminput);
            end
            
            pmin = zeros(size(LHS));
            ssmin = zeros([1 size(LHS,1)]);
            
            lb = obj.par_lb(obj.par_inds);
            ub = obj.par_ub(obj.par_inds);
            ssfun = @(pp) obj.ssfun(data,pp);
            
            for jj=1:size(LHS,1)
                p0 = LHS(jj,:);
                [pmin(jj,:) ssmin(jj)]= fmincon(ssfun,p0,[],[],[],[],lb,ub);
            end
            
            varargout{1} = pmin;
            varargout{2} = ssmin;
            
            
        end
        
        
                %%% very unfinished
        function [results, chain, s2chain, sschain] = fit_MCMC_DRAM(obj, data, initParams, options)
            
            
            
            
%             if(nargin == 0)
%                 msg = "Not enough arguments";
%                 error(msg);
%             end

            if(~exist('options','var'))
                options.nsimu = 1e4;
            end
            if(~isfield(options,'nsimu'))
                options.nsimu = 1e4;
            end
            
            if(isfield(options,'pind'))
                pind = options.pind;
                bool_pindSpecified = true;
            else
                pind = obj.par_inds;
                bool_pindSpecified = false;
            end
            
            if(isfield(options,'params_fixed'))
                validateattributes(options.params_fixed, {'struct'},{})
                params_fixed = options.params_fixed;
                if(bool_pindSpecified)
                    fixedinds = params_fixed.inds;
                else
                    fixedinds = 1:obj.par_num; %%% WRONG
                end
                params_fixed.inds;
            end
            
            
            %%% set DRAM
            
            
            %%%%% params = param_Struct_CARRGO_modified( pind, initGuess );
            params = obj.get_paramsForMCMC(initParams,pind);
            %%%%%
            
            MCMCmodel.N = sum([data.x{:}]) - 1;
            MCMCmodel.ssfun = @(p, data) obj.ssfun(data, p);
            % number of mcmc run
             if(~isfield(options,'nsimu'))
                 options.nsimu = 1*10^4;
             end

            
            ss0 = obj.ssfun(data, initParams);
            mse = ss0/(MCMCmodel.N-length(params));
            %Initial guess for error variance
            MCMCmodel.sigma2 = mse;
            
            if(~isfield(options,'updatesigma'))
                options.updatesigma = 1;  %Update s2 chain
            end

            options.stats = 1; %Include convergence statistics
            
            
            %%% run DRAM
            
            [results, chain, s2chain, sschain] = mcmcrun(MCMCmodel,data,params,options);
            

            
            
        end
            % 
        function plot_results(obj, params, data)
            
            dd.xdata = [1:.1:max(data.xdata)];
            %dd.x{1} = ~isnan( dd.xdata ); dd.x{2} = dd.x{1}; dd.ydata = zeros( size(dd.xdata,1), 2 );
            dd.x{1} = ~isnan( dd.xdata ); dd.ydata = zeros( size(dd.xdata,1), 1);
            sol_fit = obj.modelfun(dd, params);
            subplot( 2, 1, 1 ); hold on;
            
            plot( dd.xdata, sol_fit(:,1) );
            %%%% plot data
            plot( data.xdata(data.x{1}), data.ydata(data.x{1},1), 'ok', 'markerfacecolor', 'g' );
            
            
            subplot( 2, 1, 2 ); hold on;
            
            plot( dd.xdata, sol_fit(:,2) );
            %plot( dd.xdata, sol_sim(:,2), 'g' );
            %plot( data.xdata(data.x{2}), data.ydata(data.x{2},2), 'ok', 'markerfacecolor', 'g' );      
        end
        
        function multiplot_results(obj, params, data, param_names)
            
            M = size(params, 1);
            if(~exist('param_names','var'))
                param_names = strcat("Set ", string(1:M));
            end
            
            dd.xdata = [1:.1:max(data.xdata)];
            %dd.x{1} = ~isnan( dd.xdata ); dd.x{2} = dd.x{1}; dd.ydata = zeros( size(dd.xdata,1), 2 );
            dd.x{1} = ~isnan( dd.xdata ); dd.ydata = zeros( size(dd.xdata,1), 1);
            sol_fit = obj.modelfun(dd, params);
            
            
            subplot( 2, 1, 1 ); hold on;
            
            subplot( obj.y_num, 1, 1 ); hold on;
            for mm = 1:M
                sol = obj.modelfun(dd, params(mm,:));
                for yy = 1:obj.y_num
                    subplot( obj.y_num, 1, yy );
                    plot( dd.xdata, sol(:,yy) );       
                end
            end
            
            

            
           
        end

            

    end
end

