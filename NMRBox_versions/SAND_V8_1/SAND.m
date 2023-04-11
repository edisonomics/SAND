%YW wrote the base
%SAND V7: OS revised some 12-29-2022
%YW added LAMBDAMIN, SEED, LAMBDAMAX
% This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.

close all;
clear all;

%this section is only for debugging: make sure to run it under the main
%root and change parfor to for
% setenv('SANDPATH','/home/nmrbox/osanati/read-only/test1')
% setenv('ORIGIN','/home/nmrbox/osanati/read-only/usert1')
% setenv('USERF',[getenv('ORIGIN'), '/try1/'])
% setenv('USERL',[getenv('ORIGIN'), '/SAND_last_run'])
% setenv('TDD',[getenv('SANDPATH'),'/src']) 
% setenv('ELSM',[getenv('SANDPATH'),'/Edison_Lab_Shared_Metabolomics_UGA'])
% setenv('NCPUS','5')
% setenv('WDIR', getenv('USERL'))
% setenv('TEMPDIR',[getenv('USERL'),'/temp'])
% setenv('LOGDIR',[getenv('USERL'),'/log'])
% setenv('NSAMPLE','1')
% setenv('NSIG','7')  %max number of peaks in each bin'
% setenv('NTRY','2000') %number of iteration')
% setenv('TEMP','80') %temparature in MCMC
% setenv('PREP','1')
% setenv('LPPM','-1')
% setenv('UPPM','11')
% setenv('GBIN','3')
% setenv('LAMBDAMAX','15')
% setenv('LAMBDAMIN','0.01')
% setenv('SEED','1')
% mkdir /home/nmrbox/osanati/read-only/usert1/SAND_last_run/log/testveresion
% mkdir /home/nmrbox/osanati/read-only/usert1/SAND_last_run/log/testversion/results
% setenv('JOBNAME','testversion')

lastdir=getenv('USERL') %where to write last_results rundir: USERL=$ORIGIN/SAND_last_run
libdir=getenv('SANDPATH') %where scripts are
rawdir=[getenv('USERF') '/'] %where the user NMR folders data is. datadir=[paredir '/data/'];

addpath(genpath(getenv('ELSM')));
addpath(genpath(getenv('TDD')));
% addpath(genpath(getenv('WDIR')));
jobresloc2=[getenv('WDIR') '/' 'log' '/' getenv('JOBNAME') '/' 'results/']
distcomp.feature('LocalUseMpiexec',false);
tic;
shelladd=':';
defaultProfile=parallel.defaultClusterProfile;
p=parcluster(defaultProfile);
p.JobStorageLocation=getenv('TEMPDIR');
p.NumWorkers=round(str2double(getenv('NCPUS')));  
disp(p.NumWorkers)
ppool=parpool(p,p.NumWorkers - 1);
%%%%%%%%% arguments %%%%%%%%%
runid_arg=round(str2double(getenv('NSAMPLE')));
lambdamax_arg=str2double(getenv('LAMBDAMAX'));%default 15;
lambdamin_arg=str2double(getenv('LAMBDAMIN'));%default 0.01;
newprop_arg=str2double(getenv('NEWPROP'));%default 10;
thres_digit_arg=3*10^-5;
multi_replicate_arg=1;
dataset_arg=[getenv('NSAMPLE') '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runidstr=num2str(runid_arg);

dirpar=getenv('WDIR');% the last_preprocessed folder
cd([dirpar '/res/']);
system(['rm -f -r ' runidstr]);
mkdir(runidstr);
cd([runidstr]);
folder_d='deconvolution_details';
mkdir(folder_d);
datadir=[dirpar,'/res/nmrpipe_dir/',dataset_arg]


% seedi=1;
seedi=round(str2double(getenv('SEED')));
rng(seedi);
%%%%%%%%% loading data %%%%%%%%%
% spectra data for each band
locfiles_ft=strsplit(strtrim(ls([datadir 'mask_ft'])));
ft_cell=cell(1,length(locfiles_ft));
for ftfile=locfiles_ft
  [value axes]=read_nmrp([datadir 'mask_ft/' ftfile{1}]);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1;
  fttab=table(ppm,value);
  numstr=regexp(ftfile,'\d+','once','match');
  ind=str2num(numstr{1});
  ft_cell(ind)={fttab};
end
% ppm=ft_cell{1}{:,'ppm'};
car_freq_ppm=axes(8,1);%carrier frequency
resfreq=axes(14,1);%the spectrometer frequency MHZ
% mask data for each band
mask_cell=cell(1,length(locfiles_ft));
locfiles_mask=strsplit(strtrim(ls([datadir 'mask'])));
for maskfile=locfiles_mask
  [value axes]=read_nmrp([datadir 'mask/' maskfile{1}]);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1;
  masktab=table(ppm,value);
  numstr=regexp(maskfile,'\d+','once','match');
  ind=str2num(numstr{1});
  mask_cell(ind)={masktab};
end
% fid data for each band
fid_cell=cell(1,length(locfiles_ft));
locfiles_fid=strsplit(strtrim(ls([datadir 'mask_fid_matlab'])));
for fidfile=locfiles_fid
  %% load from table
  fidtab=readtable([datadir 'mask_fid_matlab/' fidfile{1}],'Format','%f%f%f');
  fidtab.Properties.VariableNames={'time','real','imag'};
  numstr=regexp(fidfile,'\d+','once','match');
  ind=str2num(numstr{1});
  fid_cell(ind)={fidtab};
end
% the ppm region table
FID=fopen([datadir 'record.txt']);
formatSpec='%s %f %f %f %f';
celltab_rang=textscan(FID,formatSpec,'MultipleDelimsAsOne',true);
fclose(FID);
rangetab=table(celltab_rang{1},celltab_rang{4},celltab_rang{5},celltab_rang{2},celltab_rang{3});%outer range, inner range
rangetab{:,1}=strrep(rangetab{:,1},'mask_fid/','');
numstr=regexp(rangetab{:,1},'\d+','once','match');
rangetab{:,1}=cellfun(@(x) {str2num(x)},numstr);
% original full spectra fid
fid_ori_tab=readtable([datadir 'ori_matlab/test_trans_ift.txt'],'Format','%f%f%f');
fid_ori_tab.Properties.VariableNames={'time','real','imag'};
fid_ori=fid_ori_tab{:,'real'}+i*fid_ori_tab{:,'imag'};
% fid before process (to calculate sigma)
fid_raw_tab=readtable([datadir 'temp.fid.txt'],'Format','%f%f%f');
fid_raw_tab.Properties.VariableNames={'time','real','imag'};
fid_raw=fid_raw_tab{:,'real'}+i*fid_raw_tab{:,'imag'};
%original full ft
[value axes]=read_nmrp([datadir 'test.ft1']);
ppm=inc2ppm(axes);
ppm=ppm.ppm1;
ft_ori_tab=table(ppm,value);
%%%%%%%%%%parameters%%%%%%%%%
% The user can reuse the following parameters. Some modifications might be needed for spectic running purpose
nregion=length(fid_cell);
region_check=1:nregion;
% sample separation
sample_ratio=[0.7 0.2 0.1];
timevec=fid_cell{1}{:,1};
ntime=length(timevec);
sample_size=floor(ntime.*sample_ratio);
sample_size(1)=ntime-sum(sample_size(2:3));
% fitting parameters
[~,cmd]=system(['csh -c "getParm -in ' datadir 'test.fid -parm FDDMXVAL"']);
shifttimeadd=round(str2num(cmd));%time point to shift. mostly 76 for experimental dataset 0 for simulated dataset
guesssig=round(str2double(getenv('NSIG')));%maximal number of possible peaks/signals
noiseseq=real(fid_raw((ntime-500):ntime));
sigmaout=str2double(getenv('SIGMA'));
if sigmaout>=0.0%use negative to represent no provding hte value
  sigma=sigmaout;
else
  sigma=std(noiseseq);%estimation of noise level
end
convfactor=1./(sample_size.*sigma^2);%for performance evaluation
nfreq=1000;%number of frequency to search
para_add_list=struct();
para_add_list.lambdadefault=15;%default value for lambda
% para_add_list.rela_range=0.1;%0.1 % to calculate searching range for parameters(when default range is not provided)
para_add_list.niteration=round(str2double(getenv('NTRY')));%number of iteration in MCMC
para_add_list.Temp=round(str2double(getenv('TEMP')));%temperature in MCMC
para_add_list.seed=seedi+1;%random seed for MCMC
para_add_list.nsig=guesssig;%maximal number of possible peaks/signals
para_add_list.train_complete=false;
para_add_list.sigma=sigma;
para_add_list.stepratio_factor=0.1;% control step size in MCMC
para_add_list.imp_ratio=0.001;% the least improvement proportion to add a new signal
para_add_list.opt_type='hybrid_freq';%optimization method
para_add_list.flag_fft=false;
para_add_list.conv_f=[car_freq_ppm resfreq];%used in conversion ppm with Hz
range_tab=[nan nan; lambdamin_arg lambdamax_arg; 0.1/10^7 20; 0 0];%f, lambda, A, phi
thres_digit=thres_digit_arg;%remove noisy regions contains no signals
para_add_list.newprop=newprop_arg;%the initial guess of new peaks need to be at least 1/newprop of the minimum among existing peaks
para_add_list.objrescale=true;%whether use sigma to rescale the objective function
% storing structures
tab_perf=zeros(nregion,3);
store_array_list=cell([length(locfiles_ft),1]);%parameter for each regions
store_spec_list=cell([length(locfiles_ft),1]);%spectra for each regions
store_fid_list=cell([length(locfiles_ft),1]);%fid for each regions
sumspec=zeros([size(ft_cell{1},1),1]);
scalespec=zeros([size(ft_cell{1},1),1]);%summed scale spectra
inforstore={};%training related information
save([folder_d '/runid' runidstr '_trainingdata.mat'],'ft_cell','mask_cell','fid_cell','rangetab','fid_ori','ft_ori_tab','para_add_list','runid_arg');%save training data
%save([jobresloc2 'runid' runidstr '_trainingdata.mat'],'ft_cell','mask_cell','fid_cell','rangetab','fid_ori','ft_ori_tab','para_add_list','runid_arg');%save training data
warning("off") %for parallel.cluster.Local
save([folder_d '/runid' runidstr '_env_beforerun.mat']);
%save([jobresloc2 'runid' runidstr '_env_beforerun.mat']);
warning("on")
%
parfor regioni=region_check%some variables are copied and modified as for parallel
  regioni
  tempinfor=struct();
  locfid=fid_cell{regioni};
  fidfile=['test' num2str(regioni) '.fid'];
  para_add_list2=para_add_list;
  para_add_list2.hdrpath=[datadir 'mask_fid/' fidfile];
  timesig=locfid{:,2}+i*locfid{:,3};
  % estimation of frequency range
  ppmrange=flip(rangetab{regioni,2:3});%the estimation range need to be larger (the whole non-zero range)
  ppmrange_sele=flip(rangetab{regioni,4:5});%the result need to be refined by inner range
  titlestr=[num2str(ppmrange(1)),'_',num2str(ppmrange(2)),'_regi',num2str(regioni)];
  frange=(ppmrange-para_add_list2.conv_f(1))*para_add_list2.conv_f(2);
  % frange_sele=(ppmrange_sele-para_add_list2.conv_f(1))*para_add_list2.conv_f(2);
  range_tab2=range_tab;
  range_tab2(1,:)=frange;
  para_add_list2.fseq=linspace(range_tab2(1,1),range_tab2(1,2),nfreq);
  para_add_list2.defaultrange=range_tab2;
  % shift for starting time
  preind=(shifttimeadd+1):length(timesig);
  fullsig=timesig;
  timesig=timesig(preind);
  timevec_sub_front=timevec(preind)-timevec(preind(1));
  % shift the fid to be cos in phase and ignore the later signals that are irregular or zeros
  [~, startime]=max(real(timesig));%start from max real
  shifttime=startime-1;
  scalespec=scalespec+mask_cell{regioni}{:,2};
  if max(diff(abs(fullsig)))<thres_digit%remove the digitial filter artifact region of noisy regions
    disp(['stop noise region' num2str(ppmrange(1)) ' ' num2str(ppmrange(2)) '_regi' num2str(regioni)]);
    continue;
  end
  resampeind=(shifttime+1):6000;%length(timesig);
  timesig=timesig(resampeind);
  timevec_sub=timevec_sub_front(resampeind);
  % multiple run and obtain the one with minimum on objective function on validaiton set
  perf_min=[];
  para_min=[];
  sampind_min=struct();
  for repi=1:multi_replicate_arg
    % the train-validtion-test separation
    % rng(repi+regioni);%seed for determinstic behaviour but will not capable of fixing behaviour for parallel
    rng(repi+regioni+para_add_list2.seed-1,'Threefry');%when need to repeat the worker results
    sampind=struct();
    allind=1:length(timesig);
    newsize=length(resampeind);
    sample_size=floor(newsize.*sample_ratio);
    sample_size(1)=newsize-sum(sample_size(2:3));
    sampind.trainind=sort(datasample(allind,sample_size(1),'Replace',false));
    sampind.validind=sort(datasample(setdiff(allind,sampind.trainind),sample_size(2),'Replace',false));
    sampind.testind=sort(setdiff(allind,[sampind.trainind,sampind.validind]));
    [para perf]=spec_est_wrap2(timesig,timevec_sub',sampind,para_add_list2);
    if repi==1
      perf_min=perf;
      para_min=para;
      sampind_min=sampind;
    else
      if perf_min.validate>perf.validate
        perf_min=perf;
        para_min=para;
        sampind_min=sampind;
      end
    end
  end
  tempinfor.sampind=sampind_min;
  tempinfor.shiftedtime=timevec_sub_front;
  perfvec=[perf_min.train perf_min.validate perf_min.test];
  perfvec=sqrt(perfvec.*convfactor);
  tab_perf(regioni,:)=perfvec;
  para_reshape=reshape(para_min,[4,floor(length(para_min)/4)])';
  [~, indsort]=sort(para_reshape(:,1));
  para_reshape=para_reshape(indsort,:);
  % para_reshape=para_reshape(para_reshape(:,1)>=frange_sele(1)&para_reshape(:,1)<=frange_sele(2),:);
  store_array_list{regioni}=para_reshape;
  % plotting each region
  newsig=sin_mixture_simu(reshape(para_min,[4,floor(length(para_min)/4)])',timevec_sub_front,nan,'complex');
  newsig=[zeros([shifttimeadd,1]); newsig];
  spec_new=ft_pipe2(table([1:length(newsig)]',real(newsig),imag(newsig)),[datadir 'test_trans.fid'],num2str(regioni),shelladd);
  new_spec_vec=spec_new{:,2};
  store_spec_list{regioni}=new_spec_vec;
  store_fid_list{regioni}=newsig;
  % fig=figure();
  % plotr(ppm,ft_cell{regioni}{:,2},'LineWidth',2);
  % hold on;
  % plotr(ppm,new_spec_vec,'LineWidth',2);%trouble might appear when the spectra is selected for regions in preprocessing; new ppm need to be used/loaded in this condition
  % set(gca,'FontSize',40);
  % xlabel('ppm');
  % xline(ppmrange_sele(1));
  % xline(ppmrange_sele(2));
  % legend('raw','estimation','bound1','bound2');
  % saveas(fig,['tempft' num2str(guesssig) '_' titlestr '_runid' runidstr '.fig'])
  % close all;
  % fig=figure();
  % plot(real(newsig));
  % hold on;
  % plot(real(fullsig));
  % saveas(fig,['tempfid' num2str(guesssig) '_' titlestr '_runid' runidstr '.fig'])
  % close all;
  % % stack spectra plot of the composed peaks
  % spec_mat_plot=[];
  % spec_mat_plot=[spec_mat_plot; ft_cell{regioni}{:,2}.'];
  % newsig=sin_mixture_simu(reshape(para_min,[4,floor(length(para_min)/4)])',timevec_sub_front,nan,'complex');
  % newsig=[zeros([shifttimeadd,1]); newsig];
  % spec_new=ft_pipe(table([1:length(newsig)]',real(newsig),imag(newsig)),[datadir 'test_trans.fid'],num2str(regioni),shelladd);
  % spec_mat_plot=[spec_mat_plot; spec_new{:,2}.'];
  % spec_mat_plot=[spec_mat_plot; spec_mat_plot(1,:)-spec_mat_plot(2,:)];
  % for sigi=1:floor(length(para_min)/4)
  %   loci=(4*(sigi-1)+1):(4*sigi)
  %   newsig=sin_mixture_simu(para_min(loci)',timevec_sub_front,nan,'complex');
  %   newsig=[zeros([shifttimeadd,1]); newsig];
  %   spec_new=ft_pipe(table([1:length(newsig)]',real(newsig),imag(newsig)),[datadir 'test_trans.fid'],num2str(regioni),shelladd);
  %   spec_mat_plot=[spec_mat_plot; spec_new{:,2}.'];
  % end
  % stackSpectra(flip(spec_mat_plot,1),ppm,0.0,10,'decompose of one nmr data set')
  % fig=gcf;
  % saveas(fig,['tempft)stack' num2str(guesssig) '_' titlestr '_runid' runidstr '.fig']);
  % close all;
  sumspec=sumspec+new_spec_vec;
  inforstore{regioni}=tempinfor;
end
save([folder_d '/runid' runidstr '_temp_store_step2.mat'],'tab_perf','store_array_list','store_spec_list','store_fid_list','sumspec','scalespec','ppm','inforstore');
%save([jobresloc2 'runid' runidstr '_temp_store_step2.mat'],'tab_perf','store_array_list','store_spec_list','store_fid_list','sumspec','scalespec','ppm','inforstore');
sumspec=sumspec./scalespec;%just to present the sum spectra which is not the final result
% visualize whole spectra (summed)
ppmrange=[-0.2,10.0];
ppmseleind=matchPPMs(ppmrange,ppm);
ppmseleseq=ppmseleind(2):ppmseleind(1);
ppmsele=ppm(ppmseleseq);
newsumsig=sumspec(ppmseleseq);
ppmbaseline=[0.1 0.3];
baserang=matchPPMs(ppmbaseline,ppmsele);
baseseq=baserang(2):baserang(1);
newsumsig=newsumsig-mean(newsumsig(baseseq));
fprintf('figure1 started')
fig=figure();
plotr(ppmsele,ft_ori_tab{ppmseleseq,2});
hold on;
plotr(ppmsele,newsumsig);
legend('FT','modeled');
saveas(fig,[folder_d '/temp_sum_comp' '_runid' runidstr '.fig']);
%saveas(fig,[jobresloc2 'temp_sum_comp' '_runid' runidstr '.fig']);
close all;

% filter and combine estimations
tabsumm=[];
% alternative region skip combination: combine results from each region and each region is separated
for regi=region_check(1:1:length(region_check))
  temptab=store_array_list{regi};
  if size(temptab,1)==0
    continue;
  end
  subreigonppm=sort([rangetab{regi,5} rangetab{regi,4}]);
  subfrange=(subreigonppm-para_add_list.conv_f(1))*para_add_list.conv_f(2);
  tabadd=temptab(temptab(:,1)>=subfrange(1)&temptab(:,1)<=subfrange(2),:);
  tabsumm=[tabsumm; tabadd];
end
preind=(shifttimeadd+1):ntime;
timevec_sub_front=timevec(preind)-timevec(preind(1));
%
paravec_propose=reshape(tabsumm',1,[]);
%
tabsumm_refine=reshape(paravec_propose,4,[])';

% % full spectra visualization
% filt_low=quantile(tabsumm(:,3),0.1);
% tabsumm_refine=tabsumm(tabsumm(:,3)>filt_low,:);
sumsig=sin_mixture_simu(tabsumm_refine,timevec_sub_front,nan,'complex');
scalfactor=0.5;
sumsig(1)=sumsig(1)*scalfactor;
sumsig=[zeros([shifttimeadd,1]); sumsig];
spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),[datadir 'test_trans.fid'],'temp',shelladd);
new_spec_vec=spec_new_sum{:,2};
fprintf('figure2 started')
fig=figure();
plotr(ppm,ft_ori_tab{:,2},'LineWidth',2);%-mean(ft_ori_tab{:,2})
hold on;

plotr(ppm,new_spec_vec,'LineWidth',2);%-mean(new_spec_vec)
legend('FT','modeled');
saveas(fig,['refine_spec_ft' '_runid' runidstr '.fig']);
%saveas(fig,[jobresloc2 'refine_spec_ft' '_runid' runidstr '.fig']);
close all;
% fid
fprintf('figure3 started')
fig=figure();
plot(real(fid_ori));
hold on;
plot(real(sumsig));
legend('raw','estimation');
saveas(fig,['refine_spec_fid' '_runid' runidstr '.fig']);
%saveas(fig,[jobresloc2 'refine_spec_fid' '_runid' runidstr '.fig']);
close all;
% objective function
obj_scaled=sqrt(sum(abs(fid_ori-sumsig).^2)/length(sumsig)/sigma^2);
% time cost
timecost=toc;
% save data
% simu_tab=table(ppm,new_spec_vec);
% simu_tab.Properties.VariableNames={'ppm','value'};
% writetable(ft_ori_tab,'raw_data.txt');
% writetable(simu_tab,'reconstruct_data.txt');
save([folder_d '/runid' runidstr '_refine_res.mat'],'timecost','obj_scaled','paravec_propose','new_spec_vec','timevec_sub_front','tabsumm_refine','para_add_list');
%save([jobresloc2 'runid' runidstr '_refine_res.mat'],'timecost','obj_scaled','paravec_propose','new_spec_vec','timevec_sub_front','tabsumm_refine','para_add_list');
warning("off") %for parallel object warning
save([folder_d '/runid' runidstr '_env_final.mat']);
%save([jobresloc2 'runid' runidstr '_env_final.mat']);
warning("on")

fnames=table2array(readtable([lastdir '/inputfolders.txt']));
for i1=1:size(fnames,1)
    if isnumeric(fnames(i1,1))
        fnamescell{i1,1}=num2str(fnames(i1,1));
    else
        fnamescell{i1,1}=fnames(i1,1);
    end

    if isnumeric(fnames(i1,2))
        fnamescell{i1,2}=num2str(fnames(i1,2));
    else
        fnamescell{i1,2}=fnames(i1,2);
    end
end
writetable(array2table(tabsumm_refine,'VariableNames',{'frequency','lambda','amplitude','phi'}),[lastdir '/res/' runidstr '/deconvoluted_peaks_refined.txt'],'Delimiter','\t');
%writetable(array2table(tabsumm_refine,'VariableNames',{'frequency','lambda','amplitude','phi'}),[jobresloc2 'deconvoluted_peaks_refined.txt'],'Delimiter','\t');

folder_r1='deconvoluted_signals_ft';
folder_r2='deconvoluted_signals_ft_txt';
folder_r3='deconvoluted_signals_fid';
folder_t='temp';
mkdir(folder_r1)
mkdir(folder_r2)
mkdir(folder_r3)
spec_new_sum=ft_pipe3(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),[datadir 'test_trans.fid'],folder_t,shelladd);
writetable(spec_new_sum,[lastdir '/res/' runidstr '/' folder_r2 '/sum_signals.txt'],'Delimiter','\t');
system(['mv ' './' folder_t '/temp.fid ' './' folder_r3 '/sum_signals.fid' ]);
system(['mv ' './' folder_t '/temp.ft ' './' folder_r1 '/sum_signals.ft' ]);
rmdir(folder_t,'s');
fig=figure();
hold on
plotr(ppm,new_spec_vec,'LineWidth',2);
for i1=1:size(tabsumm_refine,1)
    newsig=sin_mixture_simu(tabsumm_refine(i1,:),timevec_sub_front,nan,'complex');
    newsig=[zeros([shifttimeadd,1]); newsig];
    spec_new=ft_pipe3(table([1:length(newsig)]',real(newsig),imag(newsig)),[datadir 'test.fid'],folder_t,':');
    writetable(spec_new,[lastdir '/res/' runidstr '/' folder_r2 '/signal' num2str(i1) '.txt'],'Delimiter','\t');
    system(['mv ' './' folder_t '/temp.fid ' './' folder_r3 '/signal' num2str(i1) '.fid' ]);
    system(['mv ' './' folder_t '/temp.ft ' './' folder_r1 '/signal' num2str(i1) '.ft' ]);
    rmdir(folder_t,'s');
    new_spec_vec=spec_new{:,2};
    plotr(ppm,new_spec_vec,'LineWidth',1);
    
end
hold off
saveas(fig,['refined_deconvoluted' '_runid' runidstr '.fig']);
%saveas(fig,[jobresloc2 'refined_deconvoluted' '_runid' runidstr '.fig']);
close all;

for i1=1:size(fnames,1)
    if runidstr==fnamescell{i1,1}
        system(['rm -f -r ' rawdir fnamescell{i1,2} '/SAND/' getenv('JOBNAME')]);
        mkdir([rawdir fnamescell{i1,2} '/SAND/' getenv('JOBNAME')])
        system(['cp -r ' lastdir '/res/' runidstr '/* ' rawdir fnamescell{i1,2} '/SAND/' getenv('JOBNAME')]);
        mkdir([rawdir fnamescell{i1,2} '/SAND/' getenv('JOBNAME') '/preprocessed_data'])
        system(['cp -r ' lastdir '/res/nmrpipe_dir/' runidstr '/* ' rawdir fnamescell{i1,2} '/SAND/' getenv('JOBNAME') '/preprocessed_data']);
    end
end

system(['cp -r ' lastdir '/res/' runidstr '/* ' jobresloc2]);
fprintf('MATLAB run is finished')
