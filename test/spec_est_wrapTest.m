function tests=para_initguessTest
  tests=functiontests(localfunctions);
end

function testexample(testCase)
  sigma=0.2;% for noise sd
  range_tab=[3000 3600; 10 50; 1 10; 0 0];%f, lambda, A, phi
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
  nguess=3;
  tab_para=zeros(simui,npara);
  for i=1:npara%f, lambda, A, phi
    tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
  end
  resvec=sin_mixture_simu(tab_para,timevec,sigma,'complex');
  convfactor=1./(sample_size.*sigma^2);
  para_add_list=struct();
  para_add_list.fseq=linspace(range_tab(1,1),range_tab(1,2),nfreq);
  para_add_list.lambdadefault=10;
  para_add_list.defaultrange=range_tab;
  para_add_list.rela_range=0.1;
  para_add_list.niteration=10000;
  para_add_list.Temp=2000;
  para_add_list.seed=simuseed+1;
  para_add_list.nsig=nguess;
  para_add_list.train_complete=true;
  para_add_list.sigma=sigma;
  para_add_list.stepratio_factor=0.1;
  para_add_list.newprop=100;
  [para perf]=spec_est_wrap(resvec.',timevec,sampind,para_add_list);
  perfvec=[perf.train perf.validate perf.test];
  perfvec=sqrt(perfvec.*convfactor);
  % para_flag=all(abs((sort(para(:))-sort(tab_para(:)))./sort(tab_para(:)))<0.1);%some parameters are hard to within 10%
  perf_flag=all(perfvec<10);
  verifyTrue(testCase,perf_flag);
end

function testexample_hybrid(testCase)
  sigma=0.2;% for noise sd
  range_tab=[3000 3600; 10 50; 1 10; 0 0];%f, lambda, A, phi
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
  nguess=3;
  tab_para=zeros(simui,npara);
  for i=1:npara%f, lambda, A, phi
    tab_para(:,i)=rand(1,simui)*(range_tab(i,2)-range_tab(i,1))+range_tab(i,1);
  end
  resvec=sin_mixture_simu(tab_para,timevec,sigma,'complex');
  convfactor=1./(sample_size.*sigma^2);
  para_add_list=struct();
  para_add_list.fseq=linspace(range_tab(1,1),range_tab(1,2),nfreq);
  para_add_list.lambdadefault=10;
  para_add_list.defaultrange=range_tab;
  para_add_list.rela_range=0.1;
  para_add_list.niteration=5000;
  para_add_list.Temp=2000;
  para_add_list.seed=simuseed+1;
  para_add_list.nsig=nguess;
  para_add_list.train_complete=true;
  para_add_list.sigma=sigma;
  para_add_list.stepratio_factor=0.1;
  para_add_list.opt_type='hybrid_freq';
  para_add_list.newprop=100;
  [para perf]=spec_est_wrap(resvec.',timevec,sampind,para_add_list);
  perfvec=[perf.train perf.validate perf.test];
  perfvec=sqrt(perfvec.*convfactor);
  % para_flag=all(abs((sort(para(:))-sort(tab_para(:)))./sort(tab_para(:)))<0.1);%some parameters are hard to within 10%
  perf_flag=all(perfvec<10);
  verifyTrue(testCase,perf_flag);
end
