function [para perf]=spec_est_wrap(resvec,timevec,sampind,para_add_list)
% wrapper for time domain based spectra deconvolution. A iterated process will be used to eatimate parameters of underlying signal based on cross-validation.
% The program will increase the complexity (adding signals) one by one and stop when the performance didn't increase by a fixed threhold. Afterwards the model will be retrained(with local optimization) on the final model (would be the k-1 model if the k th model didn't improve performance by the threhold)
%
% Arguments:
%         resvec: vector. numeric. dimension [#time_point 1]. The time domain signal that need deconvolution. Must be provided.
%         timevec: vector. numeric. dimension [1 #time_point]. The time grid for the time domain signal. Must be provided.
%         sampind: struct. The sample index list for training, validation and testing data set. Each dimension [1 #time_point]. Must be provided.
%         para_add_list: struct. The add-on parameters for controling performance of the function. Default lambdadefault=0, rela_range=0.1
%               fseq: The frequncy sequence to sweep through. dimension [1 #f]. used in para_initial_guess. must be provided
%               lambdadefault: The default lambda value. used in para_initial_guess. default 0
%               defaultrange: The default range for all parameters. As relative step is used in the optimization, the range need to not cover 0. dimension [4 2] .default nan
%               rela_range: The relative uncerntainty for paramter searching range. used in para_initial_guess. default 0.1. This will be overwrite by defaultrange if provided. rela_range tends to be conserved when the initial estimation of A are not precise.
%               niteration: number of iterations of the MCMC. used in spec_optimize. default 1000
%               seed: random seed for MCMC. used in spec_optimize. default 1.
%               Temp: Temperature for MCMC. used in spec_optimize. default 1.
%               opt_type: string. 'default', 'hybrid', 'hybrid_freq'. Just use MCMC for optimization or use a hybrid setting (including local optimization). 'hybrid_freq' will do local optimization loc_t sweep to reduce time cost. It will also update other variables in MC steps default 'default'
%               loc_t: numeric. used only when opt_type='hybrid_freq'. Doing local optimization every loc_t sweeps. Default 100
%               stepratio_factor: step size ratio. default 0.1
%               nsig: number of expected signals. use as maximal number of signal in the sum and if -1 then no up-limit. default -1.
%               train_complete: bool. whether to continue training when performance on validation set doesn't improve with more signal. default false.
%               sigma: standard variation of noise used for early stop. default nan.
%               imp_ratio: the ratio of improvement as stop criteria. if perf_valid>(1-imp_ratio)*perf_last, then stop. default 0.1
%               flag_fft: bool. whether use fft to find the maximal A guess. If true use fft to do the initial A guess (para_initial_guess_fft), if false use (para_initial_guess). default false
%               hdrpath: string. The location of the original hdr file for fid. used when flag_fft=TRUE
%               conv_f: vector. numeric. The conversion factor from ppm to frequncy. used when flag_fft=TRUE
%               newprop: numeric. intensity of new proposed peak need to be at least min_curr_peak_intensity/para_add_list.newprop. default 10.
%               objrescale: bool. whether to use sigma to rescale the objective function. default fale
% Return:
%         para: array. numeric. The parameter for best trained model. reshape(para_min,[4,floor(length(para_min)/4)])' Each row is a different compnent and each column is a different parameter:  f, lambda, A, phi.
%         perf: struct. The performance (C_k) of the best performance model in training, validaiton and testing set.
%
% Example:
%
% n_guss_seq=1:8;%guessed number of composing sin functions
% sigma=0.2;% for noise sd
% range_tab=[3000 3600; 0.1 0.5; 1 10; -pi pi];%f, lambda, A, phi
% npara=size(range_tab,1);
% ntime=32768;
% timevec=linspace(0,1,ntime);
% simuseed=1;
% rng(simuseed);
% % cv separation
% sample_ratio=[0.8 0.1 0.1];
% sample_size=floor(ntime.*sample_ratio);
% sample_size(1)=ntime-sum(sample_size(2:3));
% sampind=struct();
% allind=1:ntime;
% sampind.trainind=datasample(allind,sample_size(1));
% sampind.validind=datasample(setdiff(allind,sampind.trainind),sample_size(2));
% sampind.testind=datasample(setdiff(allind,[sampind.trainind,sampind.validind]),sample_size(3));
% % control parameters
% nfreq=1000;
% simui=3;
% nguess=3;
% tab_para=zeros(simui,npara);
% for i=1:npara%f, lambda, A, phi
%   tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
% end
% resvec=sin_mixture_simu(tab_para,timevec,sigma,'complex');
% convfactor=1./(sample_size.*sigma^2);
% tab_perf=zeros(length(n_guss_seq),3);
% store_array=zeros(length(n_guss_seq),4*length(n_guss_seq));
% para_add_list=struct();
% para_add_list.fseq=linspace(range_tab(1,1),range_tab(1,2),nfreq);
% para_add_list.lambdadefault=0.3;
% para_add_list.defaultrange=range_tab;
% para_add_list.rela_range=0.1;
% para_add_list.niteration=10000;
% para_add_list.Temp=2000;
% para_add_list.seed=simuseed+1;
% para_add_list.nsig=nguess;
% para_add_list.train_complete=true;
% para_add_list.sigma=sigma;
% para_add_list.stepratio_factor=0.1;
% [para perf]=spec_est_wrap(resvec,timevec,sampind,para_add_list);
%
% Test:
%     results = runtests('spec_est_wrapTest.m')
% Yue Wu 06/23/2020

if ~exist('resvec','var')
  error('please provide input time domain signal');
end
if ~exist('timevec','var')
  error('please provide input time vector');
end
if ~exist('sampind','var')
  error('please provide input sample index struct for training, validaiton and testing set');
end
if ~exist('para_add_list','var')
  para_add_list=struct();
end
if ~isfield(para_add_list,'lambdadefault')
  para_add_list.lambdadefault=0;
end
if ~isfield(para_add_list,'defaultrange')
  para_add_list.defaultrange=nan;
end
if ~isfield(para_add_list,'rela_range')
  para_add_list.rela_range=0.1;
end
if ~isfield(para_add_list,'nsig')
  para_add_list.nsig=-1;
end
if ~isfield(para_add_list,'train_complete')
  para_add_list.train_complete= false;
end
if ~isfield(para_add_list,'imp_ratio')
  para_add_list.imp_ratio=0.1;
end
if ~isfield(para_add_list,'flag_fft')
  para_add_list.flag_fft=false;
end
if ~isfield(para_add_list,'newprop')
  para_add_list.newprop=10;
end
if ~isfield(para_add_list,'objrescale')
  para_add_list.objrescale=false;
end
nsig=para_add_list.nsig;

Ck=0;
para_initial_array=[];
para_iniguess_rang_comb=[];
paraarray_last=[];
para_range_last=[];
perf_last_train=[];
k=1;
resvec=resvec.';
currsig=zeros(1,length(resvec));
while true
  % parameter estimation
  disp('initial value estimation');
  if para_add_list.flag_fft
    [para_initguess para_iniguess_rang]=para_initial_guess_fft(resvec-currsig,para_add_list.hdrpath,para_add_list.conv_f,para_add_list.rela_range);
    para_initguess(3)=para_initguess(3)/length(resvec);
    % para_iniguess_rang(3,:)=para_iniguess_rang(3,:)/length(resvec);
  else
    [para_initguess para_iniguess_rang]=para_initial_guess(resvec-currsig,timevec,para_add_list.fseq,para_add_list.lambdadefault,para_add_list.rela_range);
  end
  if ~isnan(para_add_list.defaultrange)
    defaultrange=para_add_list.defaultrange;
    para_iniguess_rang=defaultrange;
    for parai=1:length(para_initguess)
      if para_initguess(parai)<defaultrange(parai,1);
        para_initguess(parai)=defaultrange(parai,1);
      end
      if para_initguess(parai)>defaultrange(parai,2);
        para_initguess(parai)=defaultrange(parai,2);
      end
    end
  end
  if k>1 & para_initguess(3)<min(para_initial_array(:,3))/para_add_list.newprop
    break;
  end
  para_initial_array=[para_initial_array; para_initguess];%k*npara, npara=4
  para_iniguess_rang_comb=[para_iniguess_rang_comb; para_iniguess_rang];%(k*npara)*2
  disp('MCMC optimization');
  [para_est perf_train]=spec_optimize(resvec(sampind.trainind),timevec(sampind.trainind),para_initial_array,para_iniguess_rang_comb,para_add_list);
  paraarray=reshape(para_est,[length(para_initguess),k])';
  para_initial_array=paraarray;
  currsig=sin_mixture_simu(paraarray,timevec,nan,'complex');
  % validation
  partial_vec=sin_mixture_simu(paraarray,timevec(sampind.validind),nan,'complex');
  perf_valid=sum(abs(resvec(sampind.validind)-partial_vec).^2);
  % if k>=2
  %   disp(['para']);
  %   disp(para_est);
  %   disp(['perf' num2str(perf_train) ' ' num2str(perf_valid)]);
  % end
  % sprintf('curr perf of %d signals %0.4f,%0.4f',k,perf_train,perf_valid)
  % disp(para_est);
  % check for improvement first and then check for maximal signal
  if (k>1)&(perf_valid>(1-para_add_list.imp_ratio)*perf_last)
    warning('adding new signal does not improve validation performance');
    if ~para_add_list.train_complete
      break;
    end
  end
  if nsig>0&k>=nsig
    warning('maximal number of signal is arrived');
    perf_last=perf_valid;
    perf_last_train=perf_train;
    paraarray_last=paraarray;
    para_range_last=para_iniguess_rang_comb;
    k=k+1;
    break;
  end
  perf_last=perf_valid;
  perf_last_train=perf_train;
  paraarray_last=paraarray;
  para_range_last=para_iniguess_rang_comb;
  k=k+1;
end
% refine training
disp('final refinement');
para_add_list_retrain=para_add_list;
para_add_list_retrain.niteration=10;
para_add_list_retrain.opt_type='hybrid';
para_add_list_retrain.Temp=1;
[para_est perf_train]=spec_optimize(resvec(sampind.trainind),timevec(sampind.trainind),paraarray_last,para_range_last,para_add_list_retrain);
paraarray=reshape(para_est,[length(para_initguess),k-1])';
% validation
partial_vec=sin_mixture_simu(paraarray,timevec(sampind.validind),nan,'complex');
perf_valid=sum(abs(resvec(sampind.validind)-partial_vec).^2);
% test
partial_vec=sin_mixture_simu(paraarray,timevec(sampind.testind),nan,'complex');
perf_test=sum(abs(resvec(sampind.testind)-partial_vec).^2);

para=para_est;
perf=struct();
perf.train=perf_train;
perf.validate=perf_valid;
perf.test=perf_test;
