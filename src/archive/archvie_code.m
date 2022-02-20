% in the following simulation, phi will be controled to constant 0
% range
wordir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/';
nsimu=5;%number of simulations
n_guss_seq=1:8;%guessed number of composing sin functions
sigma=0.2;% for noise sd
range_tab=[3000 3600; 0.1 0.5; 1 10; -7 7];%f, lambda, A, phi
npara=size(range_tab,1);
ntime=32768;
timevec=linspace(0,1,ntime);
simuseed=1;
rng(simuseed);
% cv separation
% cv=cvpartition(ones(size(timevec)),'KFold',ncv);
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
res_coll_cell={};
for simui=1:nsimu
  % simulation
  tab_para=zeros(simui,npara);
  for i=1:npara%f, lambda, A, phi
    tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
  end
  resvec=sin_mixture_simu(tab_para,timevec,sigma,'complex');
  convfactor=1./(sample_size.*sigma^2);
  tab_perf=zeros(length(n_guss_seq),3);
  store_array=zeros(length(n_guss_seq),4*length(n_guss_seq));
  for nguess=n_guss_seq
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
    paravectemp=zeros(1,4*length(n_guss_seq));
    paravectemp(1:length(para))=para;
    store_array(nguess,:)=paravectemp;
  end
  res_coll_cell=[res_coll_cell {tab_perf store_array}];
  fig=figure();
    for i=1:size(tab_perf,2)
      plot(1:nsimu,tab_perf(:,i));
    end
    legend('train','validate','test');
  saveas(fig,strcat(workdir,'deconv_examp',num2str(simui),'.fig'));
end


% single test
% simui=1;
% tab_para=zeros(simui,npara);
% for i=1:npara%f, lambda, A, phi
%   tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
% end
% resvec=sin_mixture_simu(tab_para,timevec,sigma,'real');
% f= @(x)test_objfunc(x,timevec,resvec,simui);
% lb=repmat(range_tab(1:3,1)',[1,simui]);
% ub=repmat(range_tab(1:3,2)',[1,simui]);
% otps=optimoptions('particleswarm','SwarmSize',1000);
% % [x,fval,exitflag,output]=particleswarm(f,simui*3,lb,ub,otps);
% [x]=fminsearch(f,lb);
% paratab_est=reshape(x,[3,simui]);
% paratab_est=paratab_est';
% paratab_est=[paratab_est [0 0 0]'];
% test_objfunc(x,timevec,resvec)
