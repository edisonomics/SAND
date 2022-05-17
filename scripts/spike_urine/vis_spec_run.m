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
projdir=[paredir 'result_reproduce/urine_spikin/'];
datadir=[paredir 'data/'];
libdir=[datadir 'test_trans.fid'];% a template fid file containing useful header information. https://www.dropbox.com/s/1i0dixw4vasctwu/test_trans.fid?dl=0
preresdirpath=[projdir 'res/res/deconv/res/'];
cd([projdir]);
load([projdir 'res/res/saved_preprocessing.mat']);
specppm=ppm_r;
%
sampleseq=1:nsample;
ppmrange_dss=[-0.1 0.1];
deltapm_threshold=0.002;%distance threshold for peak matching
% vis check of the deconv
for i=sampleseq
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
ft_mat_sub=[];
% simu_mat_sub=[];
% window_ppm=[-0.1 10];
% for i=1:length(delt_ppm)
%   delt=delt_ppm(i);
%   intvec=ft_mat(i,:);
%   simu_int=mat_simu_spec(i,:);
%   ppmshiftr=window_ppm-delt;%shift the ppm frame instead
%   windind=sort(matchPPMs(ppmshiftr,ppm));
%   ft_mat_ref=[ft_mat_ref; intvec(windind(1):windind(2))];
%   simu_mat_align=[simu_mat_align; simu_int(windind(1):windind(2))];
% end

windowppm=flip(window_ppm(1):((window_ppm(2)-window_ppm(1))/(size(ft_mat_ref,2)-1)):window_ppm(2));
fig=figure();
plotr(windowppm,ft_mat_ref);
savefig(fig,['ft_all_sample_ref.fig']);
close all;
fig=figure();
plotr(windowppm,simu_mat_align);
savefig(fig,['deconv_all_sample_ref.fig']);
close all;
save('deconv_res.mat','ppmlist','namelist','Alist','lamblist','binlist','cell_para','sampvec','simu_mat_align','ppm','windowppm','meta_tab_new');
save('temp_store.mat');
%
ft_mat_ref=remove_region(ft_mat_ref,ppm,4.893,4.683);
% Normalization using probabilistic quotient normalization (PQN) method
ft_mat_ref_n=normalize(ft_mat_ref,ppm,'PQN');
normcheck(ft_mat_ref_n)
fig=gcf;
saveas(fig,['normalcheck_after.fig'])
close all;
displaypeak1D(ft_mat_ref_n,ppm,0,Yvec_meta);
fig=gcf;
saveas(fig,['spectra_normalized.fig'])
close all;
ppm_XALN_PQN=vertcat(ppm,ft_mat_ref_n);
csvwrite('class_PQN.csv',ppm_XALN_PQN')
% PCA plot
%% Scale using 'logoff'
ft_mat_ref_n_s=scale(ft_mat_ref_n,'logoff');
varcheck(ft_mat_ref_n_s)


sel=[1,2];%the two groups: Control and Pure Clear Cell RCC
[Sub_XALNS2,SubT]=subset_X(XALNS2,metatab,sel,'Y');
Yvec_meta_pca=SubT.Yvec;
%% Run a Principal Component Analysis with 5 components using logoff for scaling
PCA2logoff=nipalsPCA(Sub_XALNS2,5);
close all;
VisScores(Sub_XALNS2,PCA2logoff,[1 2],'Y',Yvec_meta_pca,'conf_ellipse',false, 'showlegend', SubT.Sample_grp);
fig=gcf;
saveas(fig,['score_1_2.fig'])
close all;
VisLoadings1D(Sub_XALNS2,PCA2logoff.loadings(1,:),ppmR)
fig=gcf;
saveas(fig,['loading_1_2.fig'])
close all;

% match peaks based on ppm distance within groups
deltapm_threshold=0.02;%distance threshold for peak matching
selec_ind=find(contains(sampvec,['PD1074']));
addstr='';
%
ppmlist_selec=ppmlist(selec_ind);
namelist_selec=namelist(selec_ind);
Alist_selec=Alist(selec_ind);
lamblist_selec=lamblist(selec_ind);
binlist_selec=binlist(selec_ind);
cell_para_selec=cell_para(selec_ind);
% storage
X=ft_mat_align(selec_ind,:); % original spectral matrix
% plotr(windowppm,X);
stackSpectra(X,windowppm,-0.005,20,'alignedspec');
fig=gcf;
savefig(fig,['d1074_ref.fig']);
close all;
save(['data_to_reorder' addstr '.mat'],'X','windowppm');

% match peaks with no reorder
[ppmmatch_ind_all_nonshift ppmvec_nonshift]=ppm_list_match(ppmlist_selec_nonshift,namelist_selec_nonshift,'^unknown$',deltapm_threshold);
% update the index in the submatrix of non-reordered peaks
for sampi=1:nsample_group
  peakall=ppmlist_selec{sampi};
  mody_ind=ppmmatch_ind_all_nonshift(sampi,:);
  ind_reorder=ppmmatch_ind_all_reorder(sampi,:);
  ind_reorder=ind_reorder(~isnan(ind_reorder));
  rem_ind=sort(unique(ind_reorder));%one deconvolted peak might match tp mutiple reordered ppm by ppm distance
  lenmatch=length(rem_ind);
  reglen=diff([1 rem_ind]);
  chan_ind_vec=repelem(0:lenmatch,[reglen length(peakall)-sum(reglen)]);
  chan_ind_vec(rem_ind)=[];
  mody_ind=mody_ind+chan_ind_vec(mody_ind);
  ppmmatch_ind_all_nonshift(sampi,:)=mody_ind;
end
%
ppmmatch_ind_all=[ppmmatch_ind_all_reorder ppmmatch_ind_all_nonshift];
ppmvec=[ppmvec_reorder ppmvec_nonshift];
%
ppmmatch_ind_all_expand=ppmmatch_ind_all;
ppmvec_expand=ppmvec;
% remove NA features
naind=sum(isnan(ppmmatch_ind_all),1)>0;
ppmmatch_ind_all=ppmmatch_ind_all(:,~naind);
ppmvec=ppmvec(~naind);
%formulate other information matrix based on the ppm index
% without nan
mat_reshape_all=[];
ppm_all=[];
lambda_all=[];
cell_para_all={};
for isample=1:size(ppmmatch_ind_all,1)
  intenmat=Alist_selec{isample};
  ppmmatele=ppmlist_selec{isample};
  lammat=lamblist_selec{isample};
  locatab=cell_para_selec{isample};
  ppmind=ppmmatch_ind_all(isample,:);
  mat_reshape_all=[mat_reshape_all; intenmat(ppmind)];
  ppm_all=[ppm_all; ppmmatele(ppmind)];
  lambda_all=[lambda_all; lammat(ppmind)];
  cell_para_all=[cell_para_all; locatab(ppmind,:)];
end
% with nan
mat_reshape_all_expand=nan(size(ppmmatch_ind_all_expand));
ppm_all_expand=mat_reshape_all_expand;
for isample=1:size(ppmmatch_ind_all_expand,1)
  intenmat=Alist_selec{isample};
  ppmmatele=ppmlist_selec{isample};
  ppmind=ppmmatch_ind_all_expand(isample,:);
  no_nanind=find(~isnan(ppmind));
  mat_reshape_all_expand(isample,no_nanind)=intenmat(ppmind(no_nanind));
  ppm_all_expand(isample,no_nanind)=ppmmatele(ppmind(no_nanind));
end
% replace nan with reasonable values
for ipeak=1:size(ppmmatch_ind_all_expand,2)
  peakintvec=mat_reshape_all_expand(:,ipeak);
  nanmask=isnan(peakintvec);
  mat_reshape_all_expand(nanmask,ipeak)=min(peakintvec(~nanmask));
  ppm_all_expand(nanmask,ipeak)=mean(ppm_all_expand(~nanmask,ipeak));
end
namesall=strcat({'unknown'},cellstr(num2str(ppmvec'))');
namesall_expand=strcat({'unknown'},cellstr(num2str(ppmvec_expand'))');
save(['store_spec' addstr '.mat'],'ft_mat_align','windowppm','selec_ind','ppmmatch_ind_all','ppmvec','mat_reshape_all','ppm_all','lambda_all','cell_para_all','ppmmatch_ind_all_expand','mat_reshape_all_expand','ppm_all_expand','namesall_expand','ppmvec_expand');
save('ppm_vec.mat','ppmvec_expand','namesall_expand');
%%plotting check
samplabel=[{'ref3'},string(selec_ind')];
timevec=fid_cell{1}{:,1};
ntime=length(timevec);
shifttimeadd=76;
preind=(shifttimeadd+1):ntime;
timevec_sub_front=timevec(preind)-timevec(preind(1));
% the spectra part
spec_mat=[];
for isample=1:size(ppmmatch_ind_all_expand,1)
  tabsumm_refine=cell_para_selec{isample};
  sumsig=sin_mixture_simu(tabsumm_refine,timevec_sub_front,nan,'complex');
  scalfactor=0.5;
  sumsig(1)=sumsig(1)*scalfactor;
  sumsig=[zeros([shifttimeadd,1]); sumsig];
  spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),[datadir 'test_trans.fid'],num2str(isample));
  new_spec_vec=spec_new_sum{:,2}';
  spec_mat=[spec_mat; new_spec_vec];
end
spec_mat=[zeros([1,size(spec_mat,2)]);spec_mat];
save(['spec_recon_store' addstr '.mat'],'spec_mat');
% formulate the peak structure
% get the ppm of the show points
peakshere=struct();
for sampi=1:size(ppmmatch_ind_all_expand,1)
  spec_samp=spec_mat(1+sampi,:);
  ppm_matched=ppm_all_expand(sampi,:);
  colind=matchPPMs(ppm_matched,ppm);
  inten_matched=spec_samp(colind);
  for peaki=1:size(ppmmatch_ind_all_expand,2)
    if length(peakshere)<peaki || length(fieldnames(peakshere(peaki)))==0
      peakshere(peaki).Ridges=[-1];
      peakshere(peaki).RowInds=[1];
      peakshere(peaki).RidgeIntensities=[0];
      peakshere(peaki).CompoundNames={'unknown'};
      peakshere(peaki).quantifiable={'Y'};
    end
    peakshere(peaki).Ridges=[peakshere(peaki).Ridges ppm_matched(peaki)];
    peakshere(peaki).RowInds=[peakshere(peaki).RowInds sampi+1];
    peakshere(peaki).RidgeIntensities=[peakshere(peaki).RidgeIntensities inten_matched(peaki)];
    peakshere(peaki).CompoundNames=[peakshere(peaki).CompoundNames 'unknown'];
    peakshere(peaki).quantifiable=[peakshere(peaki).quantifiable 'Y'];
  end
end
spec_mat_flip=flip(spec_mat,1);
peakshere_flip=peakshere;
names=fieldnames(peakshere);
for peaki=1:size(peakshere,2)
  for namei=1:length(names)
    name=names{namei};
    if ~strcmp(name,'RowInds')
      peakshere_flip(peaki).(name)=flip(peakshere(peaki).(name));
    end
  end
end
samplabel_flip=flip(samplabel);
fig=stackSpectra_paintRidges_3return(spec_mat_flip,ppm,-0.00,0.005,'',peakshere_flip,10);
set(gcf,'Visible','on');
saveas(fig,['PD1074_peak_match' addstr '.fig']);
close all;
