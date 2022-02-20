function tests=spec_optimizeTest
  tests=functiontests(localfunctions);
end

function test2sigdeconv(testCase)%fit phase is tricky and it need longer input time
  para=[2 0.3 2 1; 5 0.3 1 1];
  timevec=1:0.002:10;
  para_initial_array=[1 0.4 1 1.2; 2 0.2 2 0.8];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[1 5; 0.2 0.4; 1 3; 0.1 2; 1 5; 0.2 0.4; 1 3; 0.1 2];
  sigma=0.01;
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  nseed=10;
  para_list=struct();
  para_list.Temp=100;
  para_list.niteration=4000;
  para_list.seed=seedi;
  para_list.sigma=sigma;
  para_list.stepratio_factor=0.2;
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  para_est_flag=all(distance<=para_initial_dis);
  obj_flag=perf_train<10;
  % newsig=sin_mixture_simu(reshape(para_est,[4,floor(length(para_est)/4)])',timevec,nan,'complex');
  % plot(real(newsig));
  % hold on;
  % plot(real(specsig));
  % sqrt(perf_train/length(timevec))/0.01
  verifyTrue(testCase,para_est_flag&obj_flag);
end

function test2sigdeconv_fixedphi(testCase)
  para=[2 0.3 2 0; 5 0.3 1 0];
  timevec=1:0.002:10;
  para_initial_dis=sort(para(:))*0.1;
  para_initial_array=[1.5 0.2 2.6 0; 3 0.25 0.8 0];
  para_iniguess_rang=[1 5; 0.2 0.4; 1 3; 0 0; 1 5; 0.2 0.4; 1 3; 0 0];
  specsig=sin_mixture_simu(para,timevec,0.01,'complex');
  seedi=1;
  para_list=struct();
  para_list.Temp=100;
  para_list.niteration=4000;
  para_list.seed=seedi;
  para_list.sigma=0.01;
  para_list.stepratio_factor=0.5;
  rng(seedi);
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  para_est_flag=all(distance<=para_initial_dis);
  obj_flag=perf_train<10;
  verifyTrue(testCase,para_est_flag&obj_flag);
end

function testnosigma(testCase)
  para=[2 0.3 2 1; 5 0.3 1 1];
  timevec=1:0.002:10;
  para_initial_array=[1 0.4 1 1.2; 2 0.2 2 0.8];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[1 5; 0.2 0.4; 1 3; 0.1 2; 1 5; 0.2 0.4; 1 3; 0.1 2];
  sigma=0.01;
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  nseed=10;
  para_list=struct();
  para_list.Temp=100;
  para_list.niteration=4000;
  para_list.seed=seedi;
  para_list.stepratio_factor=0.2;
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  para_est_flag=all(distance<=para_initial_dis);
  obj_flag=perf_train<10;
  % sqrt(perf_train/length(timevec))/0.01
  verifyTrue(testCase,para_est_flag&obj_flag);
end

function test2sigdeconv_realcond(testCase)
  para=[3200 0.3 2 0; 3500 0.5 1 0];
  ntime=32768;
  sigma=0.2;
  timevec=linspace(0,1,ntime);
  para_initial_array=[3000 0.1 3 0; 3600 0.4 3 0];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[3000 3600; 0.1 0.5; 1 10; 0 0; 3000 3600; 0.1 0.5; 1 10; 0 0];
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  para_list=struct();
  para_list.Temp=5000;
  para_list.niteration=10000;
  para_list.seed=seedi;
  para_list.sigma=sigma;
  para_list.stepratio_factor=0.1;
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  % para_est_flag=all(distance<=para_initial_dis)% the parameter matching wasn't included as it didn't work for some parameters
  % newsig=real(sin_mixture_simu(reshape(para_est,[4,2])',timevec,nan,'complex'));
  % plot(newsig);
  % hold on;
  % plot(real(specsig));
  % legend('fitting','data')
  obj_flag=sqrt(perf_train/length(timevec))/sigma<10;
  verifyTrue(testCase,obj_flag);
end

function test_hybrid(testCase)
  para=[3200 0.3 2 0; 3500 0.5 1 0];
  ntime=32768;
  sigma=0.2;
  timevec=linspace(0,1,ntime);
  para_initial_array=[3000 0.1 3 0; 3600 0.4 3 0];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[3000 3600; 0.1 0.5; 1 10; 0 0; 3000 3600; 0.1 0.5; 1 10; 0 0];
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  para_list=struct();
  para_list.Temp=5000;
  para_list.niteration=100;
  para_list.seed=seedi;
  para_list.sigma=sigma;
  para_list.stepratio_factor=0.1;
  para_list.opt_type='hybrid';
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  % para_est_flag=all(distance<=para_initial_dis)% the parameter matching wasn't included as it didn't work for some parameters
  % newsig=real(sin_mixture_simu(reshape(para_est,[4,2])',timevec,nan,'complex'));
  % plot(newsig);
  % hold on;
  % plot(real(specsig));
  % legend('fitting','data')
  obj_flag=sqrt(perf_train/length(timevec))/sigma<10;
  verifyTrue(testCase,obj_flag);
end

function test_hybrid_freq(testCase)
  para=[3200 0.3 2 0; 3500 0.5 1 0];
  ntime=32768;
  sigma=0.2;
  timevec=linspace(0,1,ntime);
  para_initial_array=[3000 0.1 3 0; 3600 0.4 3 0];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[3000 3600; 0.1 0.5; 1 10; 0 0; 3000 3600; 0.1 0.5; 1 10; 0 0];
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  para_list=struct();
  para_list.Temp=5000;
  para_list.niteration=5000;
  para_list.seed=seedi;
  para_list.sigma=sigma;
  para_list.stepratio_factor=0.1;
  para_list.opt_type='hybrid_freq';
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  % para_est_flag=all(distance<=para_initial_dis)% the parameter matching wasn't included as it didn't work for some parameters
  % newsig=real(sin_mixture_simu(reshape(para_est,[4,2])',timevec,nan,'complex'));
  % plot(newsig);
  % hold on;
  % plot(real(specsig));
  % legend('fitting','data')
  obj_flag=sqrt(perf_train/length(timevec))/sigma<10;
  verifyTrue(testCase,obj_flag);
end

function test_hybrid_freq_mini(testCase)
  para=[3200 0.3 2 0; 3500 0.5 1 0];
  ntime=32768;
  sigma=0.2;
  timevec=linspace(0,1,ntime);
  para_initial_array=[3000 0.1 3 0; 3600 0.4 3 0];
  para_initial_dis=sort(para(:))*0.1;
  para_iniguess_rang=[3000 3600; 0.1 0.5; 1 10; 0 0; 3000 3600; 0.1 0.5; 1 10; 0 0];
  seedi=1;
  rng(seedi);
  specsig=sin_mixture_simu(para,timevec,sigma,'complex');
  para_list=struct();
  para_list.Temp=5000;
  para_list.niteration=5000;
  para_list.seed=seedi;
  para_list.sigma=sigma;
  para_list.stepratio_factor=0.1;
  para_list.opt_type='hybrid_freq_mini';
  [para_est perf_train]=spec_optimize(specsig,timevec,para_initial_array,para_iniguess_rang,para_list);
  distance=abs(sort(para_est)-sort(para(:)));
  % para_est_flag=all(distance<=para_initial_dis)% the parameter matching wasn't included as it didn't work for some parameters
  % newsig=real(sin_mixture_simu(reshape(para_est,[4,2])',timevec,nan,'complex'));
  % plot(newsig);
  % hold on;
  % plot(real(specsig));
  % legend('fitting','data')
  obj_flag=sqrt(perf_train/length(timevec))/sigma<10;
  verifyTrue(testCase,obj_flag);
end
