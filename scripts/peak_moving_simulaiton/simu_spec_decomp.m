% an example workflow for simple simulated spectra and decompostion
% two broad peaks and one moving narrow peaks will be simulated
close all;
clear all;
%% Set your toolbox paths; functions imported from these directories:
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/Users/yuewu/Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/';
% NMR decompositon program found @ https://github.com/mikeaalv/NMR_time_domain_decomposition
localPaths.nmrdecomp_path='/Users/yuewu/Documents/GitHub/NMR_time_domain_decomposition/';
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.nmrdecomp_path));
pause(1),clc
% the path should be changed accordingly in the users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/peak_moving_simulaiton/'];
rundir=[projdir];
cd(rundir);
mkdir('res')
cd([rundir '/res/']);
seedi=1;
rng(seedi);
%%%%%%%%%%%%%%simulation parameters%%%%%%%%%%%%%%
ntime=32768;%number of time points
time1=4.6;%total time, s
freq_res=600;%observe frequency, MHZ
timevec=linspace(0,time1,ntime);
frange=[720 740 745 749 750 751 755 760];%frequency of the moving narrow peak
sigma=[0.05];%noise level
% parameter: f, lambda, amplitude, phase
range_tab=[710 790; 1 100; 0.01 10; 0 0];%range of parameters
tab_para=[760 5.0 0.05 0; 780 70.0 4 0; 750 40.0 2.2 0];%default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simuseed=1;
rng(simuseed);
%
ppm=(1:ntime)/freq_res/time1;
% plot(ppm,real(fft(resvec)));
% separation of samples
sample_ratio=[0.8 0.1 0.1];
sample_size=floor(ntime.*sample_ratio);
sample_size(1)=ntime-sum(sample_size(2:3));
sampind=struct();
allind=1:ntime;
sampind.trainind=datasample(allind,sample_size(1));
sampind.validind=datasample(setdiff(allind,sampind.trainind),sample_size(2));
sampind.testind=datasample(setdiff(allind,[sampind.trainind,sampind.validind]),sample_size(3));
% control parameters for decompostion
nfreq=1000;% the grid of frequency in initial guess% n
guesssig=4;% the guess number of signal (peaks)
nseed=3;%number of decomposition runs with different seed
para_add_list=struct();
para_add_list.fseq=linspace(range_tab(1,1),range_tab(1,2),nfreq);
para_add_list.lambdadefault=1;%default value for lambda
para_add_list.defaultrange=range_tab;
para_add_list.rela_range=0.1;% to calculate searching range for parameters(when default range is not provided)
para_add_list.niteration=2000;%number of iteration in MCMC
para_add_list.Temp=80;%temperature in MCMC
para_add_list.nsig=guesssig;%maximal number of possible peaks/signals
para_add_list.train_complete=false;
para_add_list.sigma=sigma;
para_add_list.stepratio_factor=0.1;% control step size in MCMC
para_add_list.imp_ratio=0.001;% the least improvement proportion to add a new signal
para_add_list.opt_type='hybrid_freq';%optimization method
para_add_list.newprop=500;%the initial guess of new peaks need to be at least 1/newprop of the minimum among existing peaks
%
tab_perf_record=[];
tab_para_record=[];
spec_coll=[];
parfor fragind=1:length(frange)
  rng(fragind);
  tab_para2=tab_para;
  tab_para2(1,1)=frange(fragind);
  tab_perf=[];
  store_array=[];
  temp_save_spec=[];
  % simulate combplex number fid based on the given parameters
  resvec=sin_mixture_simu(tab_para2,timevec,sigma,'complex');
  convfactor=1./(sample_size.*sigma^2);%factor for evaluating performances
  % replicates
  for seedi=1:nseed
    para_add_list2=para_add_list;
    para_add_list2.seed=seedi+1;
    % decompose the spectra
    [para perf]=spec_est_wrap(resvec.',timevec,sampind,para_add_list2);
    % store performance
    perfvec=[perf.train perf.validate perf.test];
    perfvec=sqrt(perfvec.*convfactor);
    tab_perf=[tab_perf; perfvec];
    % store parameters
    paravectemp=zeros(1,4*guesssig);
    para_reshape=reshape(para,[4,floor(length(para)/4)])';
    [~, indsort]=sort(para_reshape(:,1));
    para_reshape=para_reshape(indsort,:)';
    paravectemp(1:length(para))=para_reshape(:);
    store_array=[store_array; paravectemp];
    % visualize the estimation
    newsig=sin_mixture_simu(reshape(para,[4,floor(length(para)/4)])',timevec,nan,'complex');
    fig=figure();
    plotr(ppm,real(fft(resvec)),'LineWidth',2);
    hold on;
    plotr(ppm,real(fft(newsig)),'LineWidth',2);
    legend('simulation','estimation')
    set(gca,'FontSize',40)
    xlabel('ppm')
    saveas(fig,['widpeak_deltf' num2str(frange(fragind)) '_noise_' num2str(sigma) '_seed_' num2str(seedi) '_fft_lowalp_dense.fig']);
    close(fig);
    if seedi==1
      temp_save_spec=[temp_save_spec; [real(fft(resvec)); real(fft(newsig))]];
    end
  end
  tab_perf_record=[tab_perf_record; tab_perf];
  tab_para_record=[tab_para_record; store_array];
  spec_coll=[spec_coll; temp_save_spec];
end
save(['simu_wide_shifting_dense.mat'],'tab_perf_record','tab_para_record','spec_coll','ppm');
%plot the raw spectra
ppmrange=[1.15 1.4];
selind=1:2:15;
rangebd=sort(matchPPMs(ppmrange,ppm));
rangeseq=rangebd(1):rangebd(2);
stackSpectra(spec_coll(selind,rangeseq),ppm(rangeseq),0,100,'stack all spectra')
fig=gcf;
saveas(fig,'stackspec.fig');
close all;
% differences in f and A
fvec_real=repelem(frange,nseed);
Avec_real=repelem(0.05,size(tab_para_record,1));
[A_est,ind]=min(tab_para_record(:,[3,7,11]),[],2);%selecting all A
f_est=tab_para_record(sub2ind(size(tab_para_record),1:24,(ind'-1)*4+1));
fig=figure();
scatter(fvec_real,f_est);
saveas(fig,['f_est_scatter.fig']);
close(fig);
%
fig=figure();
scatter(Avec_real,A_est');
saveas(fig,['A_est_scatter.fig']);
close(fig);
%
fig=figure();
plotr(fvec_real/freq_res,A_est','LineWidth',3);
yline(Avec_real(1),'--','LineWidth',3);
xlabel('ppm');
ylabel('A estimation');
set(gca,'FontSize',20);
saveas(fig,['A_est_variance.fig']);
close(fig);
%
dev_perc=(Avec_real-A_est')./Avec_real;
