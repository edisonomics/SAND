function [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
% Parameter estimation for the new spectra sum for all signals and all parameters. The optimization will start from initial guess of parameter and be bound in the range. Global optimization (Metropolis Hasting) is used.
% The function can be wrong with just MCMC (default) or hybrid setting.
%           MCMC: MCMC update every variable parameters and accept/reject.
%           hybrid: MCMC update every frequency parameters, for each proposal, local optimization is implemented on f, A, and phi, the proposal is accep/reject based on objective function value on the result of local optimization. slower for the same number of iterations
%           hybrid_freq: similar to hybrid. will do local optimization loc_t sweep to reduce time cost. It will also update other variables in MC steps
%
% The matlab version of this function is slow when number of signal is large. It need to be reformulated in lower language (like C) to solve
%
% Arguments:
%         specsig: vector. numeric. The time domain signal that need deconvolution. dimension [1 #time_point]. Must be provided.
%         timevec: vector. numeric. The time grid for the time domain signal. dimension [1 #time_point]. Must be provided.
%         para_initial_array: array. numeric. The array for current estimated parameters (k-1) and initial guess of parameters (k). the initial value need not to be 0 for optimization as relative step is used.k*npara. dimension [1 4]. Must be provided.
%         para_iniguess_rang: array. numeric. The range array for all current parameters. (k*npara)*2. try not to cover 0 for this range if the corresponding parameter need optimization. dimension [4 2]. Must be provided.
%         para_list: struct. Add-on parameters.
%               niteration: number of iterations of the MCMC. used in spec_optimize. default 1000
%               seed: random seed for MCMC. used in spec_optimize. default 1.
%               Temp: Temperature for MCMC. used in spec_optimize. default 1.
%               stepratio_factor: step size ratio. default 0.1
%               sigma: standard variation of noise used for early stop. default nan.
%               opt_type: string. 'default', 'hybrid', 'hybrid_freq'. Just use MCMC for optimization or use a hybrid setting (including local optimization). 'hybrid_freq' will do local optimization loc_t sweep to reduce time cost. It will also update other variables in MC steps. 'hybrid_freq_mini' is similar to 'hybrid_freq' and have local optimization when new local optimum is arrived. default 'default'
%               loc_t: numeric. used only when opt_type='hybrid_freq'. Doing local optimization every loc_t sweeps. Default 100
%               objrescale: bool. whether to use sigma to rescale the objective function for MCMC. default fale
% Return:
%         para_est: array. numeric. All estimated paramters. Each row is a different compnent and each column is a different parameter:  f, lambda, A, phi.
%         perf_train: struct. The performance (C_k) of the model in training data set.
%
% Example:
%
% para=[2 0.3 2 0];
% timevec=1:0.002:10;
% para_initial_array=[1.5 0.2 2.6 6.28];
% para_iniguess_rang=[1 3; 0.2 0.4; 1 3; -7 7];
% specsig=sin_mixture_simu(para,timevec,0.01,'complex');
% seedi=10;
% para_list=struct();
% para_list.Temp=100;
% para_list.niteration=10000;
% para_list.seed=seedi;
% [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
%
% Test:
% results = runtests('spec_optimizeTest.m')
% Yue Wu 06/23/2020
if ~exist('specsig','var')
  error('please provide input time domain signal');
end
if ~exist('timevec','var')
  error('please provide input time grid');
end
if ~exist('para_initial_array','var')
  error('please provide input initial parameter array');
end
if ~exist('para_iniguess_rang','var')
  error('please provide parameter range array');
end
if ~exist('para_list','var')
  para_list=struct();
end
if isfield(para_list,'niteration')
  nitera=para_list.niteration;
else
  nitera=1000;
end
if isfield(para_list,'Temp')
  Temp=para_list.Temp;
else
  Temp=1;
end
if isfield(para_list,'stepratio_factor')
  stepratio_factor=para_list.stepratio_factor;
else
  stepratio_factor=0.1;
end
if isfield(para_list,'seed')
  seed=para_list.seed;
else
  seed=1;
end
if isfield(para_list,'sigma')
  sigma=para_list.sigma;
else
  sigma=nan;
end
if isfield(para_list,'opt_type')
  opt_type=para_list.opt_type;
else
  opt_type='default';
end
if isfield(para_list,'loc_t')
  loc_t=para_list.loc_t;
else
  loc_t=100;
end
if isfield(para_list,'objrescale')
  flag_scale=para_list.objrescale;
else
  flag_scale=false;
end
if flag_scale & ~isnan(sigma)
  Temp=Temp*sigma^2;
end
rng(seed);

para_initial_array_temp=para_initial_array';
paravec=para_initial_array_temp(:);
para_range=abs(para_iniguess_rang(:,2)-para_iniguess_rang(:,1));
% select out only parameters that has a nonzero range (optimizable)
para_set_opt_all=find(para_range>0);
npara_opt_all=length(para_set_opt_all);%all possible optimizable parameters
% optimizable parameters for mcmc
if strcmp(opt_type,'default')
  para_set_opt=para_set_opt_all;
  npara_opt=length(para_set_opt);
  hybriflag=0;
elseif strcmp(opt_type,'hybrid')
  para_set_opt=1:4:length(para_range);
  npara_opt=length(para_set_opt);
  hybriflag=1;
elseif strcmp(opt_type,'hybrid_freq')
  para_set_opt=para_set_opt_all;
  npara_opt=length(para_set_opt);
  hybriflag=2;
elseif strcmp(opt_type,'hybrid_freq_mini')
  para_set_opt=para_set_opt_all;
  npara_opt=length(para_set_opt);
  hybriflag=3;
end
nparadim=size(para_initial_array,2);
nsig=size(para_initial_array,1);
stepratio=stepratio_factor*para_range;

miniH=0;
paravec_min=[];
% surfcollect
vec_H=[];
array_sample_para=[];
accpe_n=0;
for sweepi=1:nitera
  % sweepi
  ind_para_ind=randi([1 npara_opt],1,npara_opt);
  ind_para=para_set_opt(ind_para_ind);
  step_sub_ratio=rand(1,length(paravec))*2-1;
  prop_gen=rand(1,npara_opt);
  for parai=1:npara_opt
    locind=ind_para(parai);
    % temppara=paravec(locind)+stepratio(locind)*step_sub_ratio(locind);%absolute step
    temppara=paravec(locind)*(1+stepratio_factor*step_sub_ratio(locind));%relative step
    if temppara<para_iniguess_rang(locind,1)
      temppara=para_iniguess_rang(locind,1);
    end
    if temppara>para_iniguess_rang(locind,2)
      temppara=para_iniguess_rang(locind,2);
    end
    paravec_propose=paravec;
    paravec_propose(locind)=temppara;
    % local optimization
    if hybriflag==1
      obj=@(x)sum_sq_obj(x,specsig,nsig,timevec,'complex');
      opts=optimoptions('fmincon','display','off');
      [paravec_propose fval]=fmincon(obj,paravec_propose,[],[],[],[],para_iniguess_rang(:,1),para_iniguess_rang(:,2),[],opts);
    end
    if hybriflag==2&mod(sweepi,loc_t)==1&mod(parai,4)==1
      obj=@(x)sum_sq_obj(x,specsig,nsig,timevec,'complex');
      opts=optimoptions('fmincon','display','off');
      [paravec_propose fval]=fmincon(obj,paravec_propose,[],[],[],[],para_iniguess_rang(:,1),para_iniguess_rang(:,2),[],opts);
    end
    partial_vec_prop=sin_mixture_simu(reshape(paravec_propose,[nparadim,nsig])',timevec,nan,'complex');
    H_prop=sum(abs(specsig-partial_vec_prop).^2);
    if sweepi==1&parai==1
      partial_vec_last=sin_mixture_simu(reshape(paravec,[nparadim,nsig])',timevec,nan,'complex');
      H_last=sum(abs(specsig-partial_vec_last).^2);
    end
    p_threh=exp(-(H_prop-H_last)/Temp);
    if prop_gen(parai)<p_threh
      paravec=paravec_propose;
      partial_vec_last=partial_vec_prop;
      H_last=H_prop;
      accpe_n=accpe_n+1;
    end
    % recording minimum values
    if sweepi==1&parai==1
      miniH=H_last;
      paravec_min=paravec;
    end
    if miniH>H_prop
      if hybriflag==3 %local exploring for new optimum
        obj=@(x)sum_sq_obj(x,specsig,nsig,timevec,'complex');
        opts=optimoptions('fmincon','display','off');
        [paravec_propose fval]=fmincon(obj,paravec_propose,[],[],[],[],para_iniguess_rang(:,1),para_iniguess_rang(:,2),[],opts);
        partial_vec_prop=sin_mixture_simu(reshape(paravec_propose,[nparadim,nsig])',timevec,nan,'complex');
        H_prop=sum(abs(specsig-partial_vec_prop).^2);
        paravec=paravec_propose;
        partial_vec_last=partial_vec_prop;
        H_last=H_prop;
      end
      miniH=H_prop;
      paravec_min=paravec_propose;
    end
    vec_H=[vec_H H_last];
  end
  array_sample_para=[array_sample_para; paravec'];
  if ~isnan(sigma)&sqrt(miniH/length(timevec))/sigma<0.1
    warning('a near optimum H was arrived')
    break;
  end
end
para_est=paravec_min;
perf_train=miniH;

% test plot on the objective function surface
% plot3(array_sample_para(:,2),array_sample_para(:,3),log(vec_H),'o');
% hist(array_sample_para(:,1));
% plot(array_sample_para(:,1),log(vec_H));
% fign=figure()
% plot(log(vec_H))
