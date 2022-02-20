% A simple single test for spec_est_wrap
wordir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/';
n_guss_seq=1:5;%guessed number of composing sin functions
sigma=0.2;% for noise sd
range_tab=[3000 3600; 0.1 0.5; 1 10; 0 0];%f, lambda, A, phi
npara=size(range_tab,1);
ntime=32768;
timevec=linspace(0,1,ntime);
simuseed=1;
rng(simuseed);
% cv separation
sample_ratio=[0.8 0.1 0.1];
sample_size=floor(ntime.*sample_ratio);
sample_size(1)=ntime-sum(sample_size(2:3));
sampind=struct();
allind=1:ntime;
sampind.trainind=datasample(allind,sample_size(1));
sampind.validind=datasample(setdiff(allind,sampind.trainind),sample_size(2));
sampind.testind=datasample(setdiff(allind,[sampind.trainind,sampind.validind]),sample_size(3));
% control parameters
nfreq=1000;
simui=3;
tab_para=zeros(simui,npara);
for i=1:npara%f, lambda, A, phi
  tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
end
resvec=sin_mixture_simu(tab_para,timevec,sigma,'complex');
convfactor=1./(sample_size.*sigma^2);
tab_perf=zeros(length(n_guss_seq),3);
store_array=zeros(length(n_guss_seq),20);
parfor nguess=n_guss_seq
  para_add_list=struct();
  para_add_list.fseq=linspace(range_tab(1,1),range_tab(1,2),nfreq);
  para_add_list.lambdadefault=0.3;
  para_add_list.lambdarangedefault=[0.1 0.5];
  para_add_list.rela_range=0.1;
  para_add_list.niteration=10000;
  para_add_list.Temp=100;
  para_add_list.seed=simuseed+1;
  para_add_list.nsig=nguess;
  para_add_list.train_complete=true;
  [para perf]=spec_est_wrap(resvec,timevec,sampind,para_add_list);
  perfvec=[perf.train perf.validate perf.test];
  perfvec=sqrt(perfvec).*convfactor;
  tab_perf(nguess,:)=perfvec;
  paravectemp=zeros(1,20);
  paravectemp(1:length(para))=para;
  store_array(nguess,:)=paravectemp;
end
save([wordir 'store_single_test_wholeworkflow.mat'],'perfvec','store_array');
fig=figure();
  for i=1:size(tab_perf,2)
    plot(1:nsimu,tab_perf(:,i));
    hold on;
  end
  legend('train','validate','test');
saveas(fig,[wordir 'single_test_wholeworkflow.fig']);
