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
% the path should be changed accordingly in hte users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/bif_urine/'];
datadir=[paredir 'data/'];
libdir=[projdir 'res/res_smallbin/proc_data/6/test_trans.fid'];% a template fid file containing useful header information. https://www.dropbox.com/s/1i0dixw4vasctwu/test_trans.fid?dl=0
preresdirpath=[projdir 'res/res_smallbin/deconv/res/'];
cd([projdir 'res/res_smallbin/']);
load([projdir 'res/res_smallbin/proc_data/saved_preprocessing.mat']);
specppm=ppm_r;
%
nsample=318;% temp
sampleseq=1:nsample;
% ppmrange_dss=[-0.1 0.1];
% deltapm_threshold=0.002;%distance threshold for peak matching
% vis check of the deconv
for i=datasample(sampleseq,10)%sampleseq
  foldpath=[preresdirpath num2str(i) '/'];
  dircont=dir([foldpath '*ft*.fig'])
  figpath=[foldpath dircont.name];
  uiopen(figpath,1);
end
close all;
obj_refine=[];
time_cost=[];
for runid=sampleseq
  load([preresdirpath num2str(runid) '/runid' num2str(runid) '_refine_res.mat']);
  obj_refine=[obj_refine obj_scaled];
  time_cost=[time_cost timecost];
end
tab_eval=table(sampleseq',obj_refine',time_cost','VariableNames',{'id','obj_refine','time_cost'});
save('performance_eval.mat','tab_eval');
% mean performance
mean(log10(tab_eval{:,'obj_refine'}))
% mean time
mean(tab_eval{:,'time_cost'})
% meta data table
readtable=readtable([projdir 'data/sample_info_all_groups.xlsx']);
% tabulate(readtable{:,'Sample_grp'})
% load raw FT of all samples
ft_mat=[];
for runid=sampleseq
  runid
  sampi_str=num2str(runid);
  load([preresdirpath sampi_str '/runid' sampi_str '_trainingdata.mat']);
  ft_mat=[ft_mat;ft_ori_tab{:,2}'];
end
fig=figure();
plotr(ppm,ft_mat);
savefig(fig,['ft_all_sample_stack.fig']);
close all;
% load deconvoluted data and formulate into list
ppmlist={};
Alist={};
lamblist={};
namelist={};
binlist={};
cell_para={};
mat_simu_spec=[];
for runid=sampleseq
  sampi_str=num2str(runid)
  load([preresdirpath sampi_str '/runid' sampi_str '_refine_res.mat']);
  load([preresdirpath sampi_str '/runid' sampi_str '_trainingdata.mat']);
  % the table: f, lambda, A, phi
  ppmvec=tabsumm_refine(:,1)/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  [ppmvec_sort sortind]=sort(ppmvec);
  ppmlist=[ppmlist; ppmvec_sort'];
  Alist=[Alist; tabsumm_refine(sortind,3)'];
  lamblist=[lamblist; tabsumm_refine(sortind,2)'];
  namelist=[namelist; {repmat({'unknown'},[1,length(ppmvec)])}];
  binlist=[binlist; rangetab{:,[4,5]}];
  cell_para=[cell_para; tabsumm_refine(sortind,:)];
  mat_simu_spec=[mat_simu_spec new_spec_vec];
end
mat_simu_spec=mat_simu_spec';
% selected interesting regions
fig=figure();
plotr(ppm,ft_mat);
savefig(fig,['ft_all_sample.fig']);
close all;
fig=figure();
plotr(ppm,mat_simu_spec);
savefig(fig,['deconv_all_sample.fig']);
close all;
save('deconv_res.mat','ppmlist','namelist','Alist','lamblist','binlist','cell_para','readtable','mat_simu_spec','ppm','ft_mat');
save('temp_store.mat','-v7.3');
%
% stack plot decompositon
exampregs=[-0.1:0.3:9.2; 0.2:0.3:9.5]';
runid=6;
esttab_loc=cell_para{runid};
ft_raw=ft_mat(runid,:);
shifttimeadd=76;
%remove noisy small signals
quantval=quantile(esttab_loc(:,3),0.2);
esttab_loc_fil=esttab_loc(esttab_loc(:,3)>quantval,:);
for regioni=1:size(exampregs,1)
  region_loc=exampregs(regioni,:);
  stackmat=[];
  % raw spectra
  stackmat=ft_raw;
  %
  ppmpara=esttab_loc_fil(:,1)/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  ppmind=find(ppmpara>region_loc(1) & ppmpara<region_loc(2));
  tabpara_loc=esttab_loc_fil(ppmind,:);
  nest=size(tabpara_loc,1);
  % sort tables
  [~,sortind]=sort(tabpara_loc(:,1));
  tabpara_loc=tabpara_loc(sortind,:);
  % the simulated spectra
  sumsig=sin_mixture_simu(tabpara_loc,timevec_sub_front,nan,'complex');
  scalfactor=0.5;
  sumsig(1)=sumsig(1)*scalfactor;
  sumsig=[zeros([shifttimeadd,1]); sumsig];
  spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,num2str(regioni));
  stackmat=[stackmat; spec_new_sum{:,2}'];
  % each deconv component
  for paraseti=1:nest
    sumsig=sin_mixture_simu(tabpara_loc(paraseti,:),timevec_sub_front,nan,'complex');
    scalfactor=0.5;
    sumsig(1)=sumsig(1)*scalfactor;
    sumsig=[zeros([shifttimeadd,1]); sumsig];
    spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,num2str(paraseti));
    stackmat=[stackmat; spec_new_sum{:,2}'];
  end
  % color settings
  colorset=struct();
  colorset.rgb=flip([[0 0 0]; [1 0 0]; repmat([0 0 0.7],[nest,1])],1);
  colorset.categories=table(flip([{'ft'},{'sum'},{'estimation'}]'));
  colorset.colorList=flip([0 0 0; 1 0 0; 0 0 0.7],1);
  %
  ppmvis_rang=sort(matchPPMs(region_loc,ppm));
  ppmvis_ind=ppmvis_rang(1):ppmvis_rang(2);
  stackmat=flip(stackmat,1);
  stackSpectra(stackmat(:,ppmvis_ind),ppm(ppmvis_ind),0.0,1,['example ' num2str(regioni)],'colors',colorset);
  fig=gcf;
  saveas(fig,['stack_example_region_' num2str(regioni) '.fig']);
  close all;
end
