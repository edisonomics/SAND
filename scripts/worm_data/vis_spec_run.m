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
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/result_reproduce/worm_data/'
projdir=[paredir 'res/res/'];
cd(projdir);
dirdeconv=[projdir 'deconv/res/'];
datadir=[paredir '/data/'];
libdir=[datadir 'test_trans.fid']
nsample=146;
meta_tab=readtable([projdir 'metatab_all.txt'],'Delimiter','\t');
sampvec=meta_tab{:,'genotype'};
% load example ft
sampi_str='3';
load([dirdeconv sampi_str '/runid' sampi_str '_trainingdata.mat']);
ft_examp=ft_ori_tab;
ppm=ft_examp{:,1};
%
sampseq=[1:21 23:146];%1:nsample; one sample omitted in the analysis
sampvec=sampvec(sampseq);
meta_tab_new=meta_tab(sampseq,:);
% load raw FT of all samples
ft_mat=[];
for sampi=sampseq
  sampi
  sampi_str=num2str(sampi);
  load([dirdeconv sampi_str '/runid' sampi_str '_trainingdata.mat']);
  ft_mat=[ft_mat;ft_ori_tab{:,2}'];
end
fig=figure();
plotr(ppm,ft_mat);
savefig(fig,['ft_all_sample_stack.fig']);
close all;
% ref ppm dss range
rangedss=[-0.1 0.25];
dssind=sort(matchPPMs(rangedss,ppm));
[intmax maxind]=max(ft_mat(:,dssind(1):dssind(2)),[],2);
delt_ppm=-ppm(dssind(1)+maxind-1);%shift for spectra
% load deconvoluted data and formulate into list
ppmlist={};
Alist={};
lamblist={};
namelist={};
binlist={};
cell_para={};
mat_simu_spec=[];
for sampi_i=1:length(sampseq)
  sampi=sampseq(sampi_i)
  sampi_str=num2str(sampi);
  load([dirdeconv sampi_str '/runid' sampi_str '_refine_res.mat']);
  load([dirdeconv sampi_str '/runid' sampi_str '_trainingdata.mat']);
  % shift f
  ppmtemp=tabsumm_refine(:,1)/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  ppmtemp=ppmtemp+delt_ppm(sampi_i);
  tabsumm_refine(:,1)=(ppmtemp-para_add_list.conv_f(1))*para_add_list.conv_f(2);
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
% reference original spectra by DSS
ft_mat_align=[];
simu_mat_align=[];
window_ppm=[-0.1 10];
for i=1:length(delt_ppm)
  delt=delt_ppm(i);
  intvec=ft_mat(i,:);
  simu_int=mat_simu_spec(i,:);
  ppmshiftr=window_ppm-delt;%shift the ppm frame instead
  windind=sort(matchPPMs(ppmshiftr,ppm));
  ft_mat_align=[ft_mat_align; intvec(windind(1):windind(2))];
  simu_mat_align=[simu_mat_align; simu_int(windind(1):windind(2))];
end
windowppm=flip(window_ppm(1):((window_ppm(2)-window_ppm(1))/(size(ft_mat_align,2)-1)):window_ppm(2));
fig=figure();
plotr(windowppm,ft_mat_align);
savefig(fig,['ft_all_sample_ref.fig']);
close all;
fig=figure();
plotr(windowppm,simu_mat_align);
savefig(fig,['deconv_all_sample_ref.fig']);
close all;
save('deconv_res.mat','ppmlist','namelist','Alist','lamblist','binlist','cell_para','sampvec','simu_mat_align','ppm','windowppm','meta_tab_new');
save('temp_store.mat');
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
savefig(fig,['ft_pd1074_ref.fig']);
close all;
save(['data_to_reorder' addstr '.mat'],'X','windowppm');

% spectral reordering based 'alignment'
load([datadir 'reorder_output.mat']);
ppmmat=AlignedPeaks_o.shift+AlignedPeaks_o.shiftV;
[~,sortind]=sort(mean(ppmmat,1));
ppmmat=ppmmat(:,sortind);
% match peaks with the peak reorder matrix
norderpeak=size(ppmmat,2);
nsample_group=length(ppmlist_selec);
ppmlist_selec_nonshift=ppmlist_selec;
namelist_selec_nonshift=namelist_selec;
%
ppmvec_reorder=[];
ppmmatch_ind_all_reorder=[];
flagfini=ones([norderpeak,1]);
% match deconv ppm to reorder ppm
for peaki=1:norderpeak
  ppm_order=ppmmat(:,peaki);
  ppmmatchind=[];
  % find corresponding match in all samples
  ppmcoll=[];
  for sampi=1:nsample_group
    ppmvec=ppmlist_selec_nonshift{sampi};
    [ppmdist minind]=min(abs(ppmvec-ppm_order(sampi)));
    if ppmdist<deltapm_threshold/2
      ppmmatchind=[ppmmatchind minind];
      ppmcoll=[ppmcoll ppmvec(minind)];
    else %some reordered peaks might not match to deconvoluted peaks
      ppmmatchind=[ppmmatchind NaN];
      flagfini(peaki)=0;
      % warning('missing peak matching between reorder and deconvolution');
    end
  end
  ppmvec_reorder=[ppmvec_reorder mean(ppmcoll)];
  ppmmatch_ind_all_reorder=[ppmmatch_ind_all_reorder ppmmatchind'];
end
% ppm list for non-reorder peaks
% update for complete matches (remove incomplete ones)
for sampi=1:nsample_group
  matchind=ppmmatch_ind_all_reorder(sampi,:);
  matchind=matchind(~isnan(matchind));
  ppmlist_selec_nonshift{sampi}(matchind)=[];
  namelist_selec_nonshift{sampi}(matchind)=[];
end
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
%% correlation based network
threhold_corr=0.9;
cor_intensity=corr(mat_reshape_all_expand,mat_reshape_all_expand,'Type','Pearson');
indmat=zeros(size(cor_intensity));
lenmat=size(cor_intensity,1);
for i=1:lenmat
  indmat(i,(i+1):lenmat)=1;
end
triaind=find(indmat);
cor_intensity_triag_vec=cor_intensity(triaind);
cor_intensity_triag=cor_intensity;
cor_intensity_triag(~indmat)=0;
%%edge threhold
thres_display_vec_inten=[quantile(cor_intensity_triag_vec,threhold_corr)];
%%%positive correlation of intensity
[rowind colind]=find(cor_intensity_triag>=thres_display_vec_inten);
node1={};
node2={};
correlation=[];
for elei=1:length(rowind)
  roweleind=rowind(elei);
  coleleind=colind(elei);
  node1{elei}=namesall_expand{roweleind};
  node2{elei}=namesall_expand{coleleind};
  correlation(elei)=cor_intensity_triag(roweleind,coleleind);
end
cornet_table_corr=table(node1',node2',correlation','VariableNames',{'source' 'target' 'association'});
writetable(cornet_table_corr,['deconv_corr_PD1074.0_9.pearson.norm' addstr '.txt'],'Delimiter','\t');
% load the correlation network in cytoscape, visualize, mcl cluster with granulity 4, output the node table
%plot clusters of correlation network with spectra
clutstr='0_9';
tabclust=readtable([projdir 'deconv_corr_PD1074.0_9.pearson.norm.txt--clustered_mcl4node.csv']);
clustervec=tabclust{:,'x__mclCluster'};
% single nodes will be ignored by tabulate
clust_count=cell2table(tabulate(clustervec),'VariableNames',{'Value','Count','Percent'});
clustertab_sele=clust_count(clust_count{:,'Count'}>=3,:);
clustertab_sele=sortrows(clustertab_sele,'Count','descend');
cluster_sele=clustertab_sele{:,'Value'};
%
plotstore_cell={};
save('temp_store2.mat');
for clusti=1:length(cluster_sele)
  clusthere=cluster_sele{clusti};
  nodenames=tabclust{strcmp(tabclust{:,'x__mclCluster'},clusthere),'name'};
  nodeind=cellfun(@(x) find(strcmp(namesall_expand,x)),nodenames);
  recordind=ppmmatch_ind_all_expand(:,nodeind);
  spec_mat=[];
  for sampi=1:size(ppmmatch_ind_all_expand,1)
    paratab=cell_para_selec{sampi};
    loctabind=recordind(sampi,:);
    loctabind=loctabind(~isnan(loctabind));
    if length(loctabind)==0
      continue;
    end
    loctab=paratab(loctabind,:);
    sumsig=sin_mixture_simu(loctab,timevec_sub_front,nan,'complex');
    scalfactor=0.5;
    sumsig(1)=sumsig(1)*scalfactor;
    sumsig=[zeros([shifttimeadd,1]); sumsig];
    spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),[datadir 'test_trans.fid'],num2str(sampi));
    new_spec_vec=spec_new_sum{:,2}';
    spec_mat=[spec_mat; new_spec_vec];
  end
  close all;
  stackSpectra(spec_mat,ppm,0.0,10,['cluster set ' clusthere]);
  fig=gcf;
  saveas(fig,['stack_cluster_' num2str(clusthere) '.' clutstr '.fig']);
  close all;
  plotstore_cell{clusti}=spec_mat;
end
save(['cor_clut_spec_plotdata.' clutstr '.mat'],'plotstore_cell');

% visualize some clusters with matched compounds
ppmloc=ppm;
% load GISSMO spectra
%
compounds=[{'R-Lactate'},{'L-Isoleucine'}];
clusters=[{{'9'}},{{'17'}}];
ppmrange=[{[{[1.25 1.4]},{[4 4.2]}]},{{[0.9 1.05]}}];
cluster_match=table(compounds',clusters',ppmrange','VariableNames',{'compounds' 'clusters' 'ppmrange'});
stackmat=[];
% the whole spectra
runid=1;
esttab_loc=cell_para_selec{runid};
sumsig=sin_mixture_simu(esttab_loc,timevec_sub_front,nan,'complex');
scalfactor=0.5;
sumsig(1)=sumsig(1)*scalfactor;
sumsig=[zeros([shifttimeadd,1]); sumsig];
spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,'temp');
stackmat=[stackmat; spec_new_sum{:,2}'];
lensize=[];
for rowi=1:size(cluster_match,1)
  name=cluster_match{rowi,'compounds'};
  clusters=cluster_match{rowi,'clusters'}{1};
  region_loc=cluster_match{rowi,'ppmrange'}{1};
  indtab=find(strcmp(name,tabinfor{:,'CompoundName'}));
  id=tabinfor{indtab,'EntryID'};
  match_spec_id=[id{1} '_simulation_1'];%use the information from the first simulation in gissmo for each compound
  matchspec=strdata.(match_spec_id);%compound spectra
  %
  % cluster peaks
  tempmat=[];
  for clusthere=clusters
    nodenames=tabclust{strcmp(tabclust{:,'x__mclCluster'},clusthere),'name'};
    nodeind=cellfun(@(x) find(strcmp(namesall_expand,x)),nodenames);
    recordind=ppmmatch_ind_all_expand(runid,nodeind);
    recordind=recordind(~isnan(recordind));
    if length(recordind)==0
      continue;
    end
    loctab=esttab_loc(recordind,:);
    sumsig=sin_mixture_simu(loctab,timevec_sub_front,nan,'complex');
    scalfactor=0.5;
    sumsig(1)=sumsig(1)*scalfactor;
    sumsig=[zeros([shifttimeadd,1]); sumsig];
    spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,'temp');
    tempmat=[tempmat; spec_new_sum{:,2}'];
  end
  tempmat_sum=sum(tempmat,1);
  %
  colorset=struct();
  colorset.rgb=flip([[0 0 0]; [0 0 1]; [0 1 0]],1);
  colorset.categories=table(flip([{'sum'},{'cluster'},{'reference'}]'));
  colorset.colorList=flip([0 0 0; 0 0 1; 0 1 0;],1);
  % plot regions
  for regionhere=region_loc
    regionhere=regionhere{1};
    ppmvis_rang=sort(matchPPMs(regionhere,ppmloc));
    ppmvis_ind=ppmvis_rang(1):ppmvis_rang(2);
    %
    matchspec_inte=interp1(ppm',matchspec,ppmloc');
    matchspec_inte=matchspec_inte./max(matchspec_inte(ppmvis_ind))*max(tempmat_sum(ppmvis_ind));
    stackmat2=stackmat;
    stackmat2=[stackmat2; tempmat_sum; matchspec_inte];
    stackmat2=flip(stackmat2,1);
    stackSpectra(stackmat2(:,ppmvis_ind),ppmloc(ppmvis_ind),0.0,10,[name{1}],'colors',colorset);
    fig=gcf;
    saveas(fig,['stack_clusters_' name{1} num2str(regionhere(1)) '_' num2str(regionhere(2)) '.fig']);
    close all;
  end
  lensize=[lensize size(tempmat_sum,1)+1];
end
