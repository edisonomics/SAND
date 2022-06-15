close all;
clear all;
%% Set your toolbox paths; functions imported from these directories:
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/Users/yuewu/Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/';
% NMR decompositon program found @ https://github.com/edisonomics/SAND
localPaths.nmrdecomp_path='/Users/yuewu/Documents/GitHub/SAND/';
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.nmrdecomp_path));
pause(1),clc
% the path should be changed accordingly in hte users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/urine_like_simulaiton_group/'];
datadir=[paredir 'data/urine_fitting/'];
libdir=[datadir 'test_trans.fid'];% a template fid file containing useful header information. https://www.dropbox.com/s/1i0dixw4vasctwu/test_trans.fid?dl=0
preresdirpath=[projdir 'res/deconv/res/res/'];
cd([projdir]);
% load original simulation information
load([projdir 'res/saved_simulation.mat']);
specppm=ppm_r;
%
sampleseq=1:nsample;
ppmrange_dss=[-0.1 0.1];
deltapm_threshold=0.002;%distance threshold for peak matching
%
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

% remove broad peaks
groundtruth_tab=groundtruth_tab(groundtruth_tab{:,'lambda'}<=15,:);
% ground truth
for sampi=sampleseq
  rowinds=find(groundtruth_tab{:,'simulation'}==sampi);
  temptab=groundtruth_tab(rowinds,:);
  temptab{:,'frequency'}=temptab{:,'frequency'}/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  regind=temptab{:,1}>ppmrange_dss(1) & temptab{:,1}<ppmrange_dss(2);
  dssconc=max(temptab{regind,'A'});
  temptab{:,'A'}=temptab{:,'A'}/dssconc;%normalize to DSS
  groundtruth_tab(rowinds,:)=temptab;
end
groundtruth_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation','group'};

% load decomposation estimation of different spectra
namelist={};
est_tab=[];%PPM, lambda, A, simulation_ind
for sampi=1:nsample
  sample=sampleseq(sampi);
  samplestr=num2str(sample);
  load([preresdirpath samplestr '/runid' samplestr '_env_final.mat']);
  runtab=array2table(tabsumm_refine,'VariableNames',{'PPM','lambda','A','phase'});%f, lambda, A, phi
  nfeature=size(runtab,1);
  runtab{:,'PPM'}=runtab{:,'PPM'}/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  % DSS intenstiy
  regind=runtab{:,'PPM'}>ppmrange_dss(1) & runtab{:,'PPM'}<ppmrange_dss(2);
  dssconc=max(runtab{regind,'A'});
  runtab{:,'A'}=runtab{:,'A'}/dssconc;%normalize to DSS
  simu=repmat(sampi,[nfeature,1]);
  runtab=[runtab table(simu)];
  namelist=[namelist; {repmat({'unknown'},[1,nfeature])}];
  est_tab=[est_tab; runtab];
end
est_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};

%intensity and integral based estimation
est_other_tab=[];
for sampi=1:size(specmat,1)
  spec_here=specmat(sampi,:);
  runtab=[];
  for bini=1:size(binrange,1)
    ppmrange=binrange(bini,:);
    [indrang]=sort(matchPPMs(ppmrange,specppm));
    indseq=indrang(1):indrang(2);
    spec_reg_shift=spec_here(indseq);
    baseval=min(spec_reg_shift);
    % baseval=0;
    [est_inten maxind]=max(spec_reg_shift-baseval);
    est_auc=trapz(spec_reg_shift-baseval);
    est_temp_tab=table(specppm(indrang(1)+maxind-1),nan,est_inten,est_auc,0,sampi);
    runtab=[runtab; est_temp_tab];
  end
  regind=runtab{:,1}>ppmrange_dss(1) & runtab{:,1}<ppmrange_dss(2);
  %normalize to DSS
  runtab{:,3}=runtab{:,3}/max(runtab{regind,3});
  runtab{:,4}=runtab{:,4}/max(runtab{regind,4});
  est_other_tab=[est_other_tab; runtab];
end
est_other_tab.Properties.VariableNames={'PPM','lambda','intensity','integral','phase','simulation'};
%
quan_str=struct();
quan_str.deconv=est_tab;
quan_str.intensity=est_other_tab(:,{'PPM','lambda','intensity','phase','simulation'});
quan_str.integral=est_other_tab(:,{'PPM','lambda','integral','phase','simulation'});
% select within interested ppm range
selerange=[6.8 8.6; -0.1 0.1];% considered ppm range
for type=fieldnames(quan_str)'
  type=type{1};
  tempdata=quan_str.(type);
  rowind_ppm=find(tempdata{:,'PPM'}>selerange(1,1)&tempdata{:,'PPM'}<selerange(1,2) | tempdata{:,'PPM'}>selerange(2,1)&tempdata{:,'PPM'}<selerange(2,2));
  tempdata=tempdata(rowind_ppm,:);
  quan_str.(type)=tempdata;
end
% match estimation with ground truth
summ_str=struct();
for type=fieldnames(quan_str)'
  type=type{1};
  summtab=[];
  rec_ratio=[];
  est_tab_temp=quan_str.(type);
  est_tab_temp.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};
  for simui=sampleseq
    subtab_est=est_tab_temp(est_tab_temp{:,'simulation'}==simui,:);
    subtab_true=groundtruth_tab(groundtruth_tab{:,'simulation'}==simui,:);
    distmat=abs(pdist2(subtab_est{:,'PPM'},subtab_true{:,'PPM'}));
    %%%selecting out ppm points that are pairwise closest to each other
    [ppm_match_val1,ppm_match_ind1]=min(distmat,[],2);
    [ppm_match_val2,ppm_match_ind2]=min(distmat,[],1);
    ind_true=[];
    ind_est=[];
    ppm_match_val=[];
    for ppm_min_i=1:length(ppm_match_ind1)
      if ppm_match_ind2(ppm_match_ind1(ppm_min_i))==ppm_min_i%check for pairwise match
        ppm_match_val=[ppm_match_val ppm_match_val1(ppm_min_i)];
        ind_true=[ind_true ppm_match_ind1(ppm_min_i)];
        ind_est=[ind_est ppm_min_i];
      end
    end
    % filter by ppm distances
    thres_ind=find(ppm_match_val<deltapm_threshold);
    ind_est=ind_est(thres_ind);
    ind_true=ind_true(thres_ind);
    subtab_est_match=subtab_est(ind_est,{'PPM','A','lambda','phase'});
    subtab_est_match.Properties.VariableNames={'PPM_est','A_est','lambda_est','phase_est'};
    subtab_true_match=subtab_true(ind_true,{'PPM','A','lambda','phase','group','simulation'});
    subtab_true_match.Properties.VariableNames={'PPM_true','A_true','lambda_true','phase_true','group','simulation'};
    summtab=[summtab; [subtab_est_match subtab_true_match]];
    rec_ratio=[rec_ratio size(subtab_est_match,1)/size(subtab_true,1)];
  end
  % summtab=summtab(summtab{:,'A_est'}>10^-3 & summtab{:,'A_true'}>10^-3,:);
  tempstr=struct();
  tempstr.summtab=summtab;
  tempstr.rec_ratio=rec_ratio;
  summ_str.(type)=tempstr;
end
% calculate evaluations
evalu_str=struct();
for type=fieldnames(quan_str)'
  type=type{1};
  summtab=summ_str.(type).summtab;
  rel_mse_vec=[];
  mse_vec=[];
  corxy_vec=[];
  k_vec=[];
  for simui=sampleseq
    loctab=summtab(summtab{:,'simulation'}==simui,:);
    xvec=loctab{:,'A_true'};
    yvec=loctab{:,'A_est'};
    ndata=length(xvec);
    rel_mse_vec=[rel_mse_vec sum(((xvec-yvec)./mean([xvec yvec],2)).^2)/ndata];
    mse_vec=[mse_vec sum((xvec-yvec).^2)/ndata];
    corxy_vec=[corxy_vec corr(xvec,yvec)];
    dlm=fitlm(xvec,yvec,'Intercept',false);
    k_vec=[k_vec dlm.Coefficients.Estimate];
  end
  evalu=struct();
  for eval_ele={'rel_mse' 'mse' 'corxy' 'k'}
    eval_ele=eval_ele{1};
    locvec=eval([eval_ele '_vec']);
    evalu.(eval_ele)=mean(locvec);
    evalu.([eval_ele '_ste'])=std(locvec)/sqrt(length(locvec));
  end
  evalu_str.(type)=evalu;
end
% make the table
evalu_tab=cell2table(cell(0,5),'VariableNames',{'rel_mse','mse','corxy','k','quan_method'});
for methele=fieldnames(evalu_str)'
  methele=methele{1};
  loctab=struct2table(evalu_str.(methele));
  loctab=[loctab(:,{'rel_mse','mse','corxy','k'}) table({methele},'VariableNames',{'quan_method'})];
  evalu_tab=[evalu_tab; loctab];
end
save('evaluation.mat','evalu_str','evalu_tab');
% scattter plot
for type=fieldnames(summ_str)'
  type=type{1};
  summtab=summ_str.(type).summtab;
  evalu=evalu_str.(type);
  h=figure();
    gscatter(summtab{:,'A_true'},summtab{:,'A_est'},summtab{:,'simulation'},[],[],[20]);
    xlabel('ground truth');
    ylabel('estimation');
    title([type ' correlation ' num2str(evalu.corxy), ' mse ' num2str(evalu.mse) ' k ' num2str(evalu.k)]);
  saveas(h,['scatter_simulation.' type '.fig']);
  close(h);
end
% peaks that are in simulation but not recoverd.
rec_ratio_vec=[];
for type=fieldnames(summ_str)'
  type=type{1};
  rec_ratio_vec=[rec_ratio_vec mean(summ_str.(type).rec_ratio)];
end
% lambda estimation
summtab=summ_str.deconv.summtab;
h=figure();
  gscatter(summtab{:,'lambda_true'},summtab{:,'lambda_est'},summtab{:,'simulation'},[],[],[20]);
  xlabel('ground truth');
  ylabel('estimation');
  title([' lambda ']);
saveas(h,['scatter_simulation.' '_lambda.fig']);
close(h);
corr(summtab{:,'lambda_true'},summtab{:,'lambda_est'})

% correlation network
%% match peaks from different samples
ppmlist={};
namelist={};
intlist={};
for runid=sampleseq
  loctab=est_tab(est_tab{:,'simulation'}==runid,:);
  ppmlist=[ppmlist loctab{:,'PPM'}'];
  namelist=[namelist {repmat({'unknown'},[1,size(loctab,1)])}];
  intlist=[intlist loctab{:,'A'}];
end
[ppmmatch_ind_all ppmvec]=ppm_list_match(ppmlist',namelist,'^unknown$',deltapm_threshold);
mat_reshape_all=[];
for isample=1:size(ppmmatch_ind_all,1)
  intenvec=intlist{isample};
  ppmind=ppmmatch_ind_all(isample,:);
  intenarray=intenvec(ppmind);
  mat_reshape_all=[mat_reshape_all; intenarray'];
end
%% match with ground truth peak sets
subtab_true=groundtruth_tab(groundtruth_tab{:,'simulation'}==1,:);%all simulation has the same peak set
distmat=abs(pdist2(ppmvec',subtab_true{:,'PPM'}));
[ppm_match_val1,ppm_match_ind1]=min(distmat,[],2);
[ppm_match_val2,ppm_match_ind2]=min(distmat,[],1);
ind_true=[];
ind_est=[];
for ppm_min_i=1:length(ppm_match_ind1)
  if ppm_match_ind2(ppm_match_ind1(ppm_min_i))==ppm_min_i%check for pairwise match
    ind_true=[ind_true ppm_match_ind1(ppm_min_i)];
    ind_est=[ind_est ppm_min_i];
  end
end
mat_reshape_all=mat_reshape_all(:,ind_est);
ppmvec=ppmvec(ind_est);
% remove DSS as it cause numeric problem as its not changing
remind=find(abs(ppmvec)<0.01);
ppmvec(remind)=[];
mat_reshape_all(:,remind)=[];
%
namesall=strcat({'unknown'},cellstr(num2str(ppmvec'))');
%% correlation network
threhold_corr=0.9;
cor_intensity=corr(mat_reshape_all,mat_reshape_all,'Type','Spearman');
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
thres_display_vec_inten=0.9
%%%positive correlation of intensity
[rowind colind]=find(cor_intensity_triag>=thres_display_vec_inten);
node1={};
node2={};
correlation=[];
for elei=1:length(rowind)
  roweleind=rowind(elei);
  coleleind=colind(elei);
  node1{elei}=namesall{roweleind};
  node2{elei}=namesall{coleleind};
  correlation(elei)=cor_intensity_triag(roweleind,coleleind);
end
cornet_table_corr=table(node1',node2',correlation','VariableNames',{'source' 'target' 'association'});
writetable(cornet_table_corr,'deconv_corr_thres0_9.txt','Delimiter','\t');
% load and visualize in cytoscape
% find connected groups from correaltion edges
% if the cluster pattern is already clear without community clustering
cormat_thres=cor_intensity>=thres_display_vec_inten;
G=graph(cormat_thres,'omitselfloops');
% plot(G);
[graph_conn_group]=conncomp(G);
groudtruth_onesamp=groundtruth_tab(groundtruth_tab{:,'simulation'}==1,:);
%
ppm_true=groudtruth_onesamp{:,'PPM'};
distmat=abs(pdist2(ppmvec',ppm_true));
[ppm_match_val1,ppm_match_ind1]=min(distmat,[],2);
[ppm_match_val2,ppm_match_ind2]=min(distmat,[],1);
ind_true=[];
ind_est=[];
for ppm_min_i=1:length(ppm_match_ind1)
  if ppm_match_ind2(ppm_match_ind1(ppm_min_i))==ppm_min_i%check for pairwise match
    ind_true=[ind_true ppm_match_ind1(ppm_min_i)];
    ind_est=[ind_est ppm_min_i];
  end
end
%
matchvec=zeros(size(namesall));
matchratio=[];
matchcluster=[];
for groupi=1:(ncompound-1)%ignore DSS
  groupind=find(groudtruth_onesamp{:,'group'}==groupi);
  groupsize=length(groupind);
  matchedind=find(ismember(ind_true,groupind));
  coll_clust=[];
  matchcoll={};
  for estclut=unique(graph_conn_group)
    groupind=find(graph_conn_group==estclut);
    match_est_ind=ind_est(matchedind);
    memmatch=find(ismember(match_est_ind,groupind));
    coll_clust=[coll_clust length(memmatch)];
    matchcoll=[matchcoll {match_est_ind(memmatch)}];
  end
  [maxsize maxind]=max(coll_clust);
  matchratio=[matchratio maxsize/groupsize];
  matchvec(matchcoll{maxind})=1;
  matchcluster=[matchcluster maxind];
end
fig=figure();
histogram(matchratio);
xlabel('cluster recover ratio');
ylabel('frequency');
set(gca,'FontSize',40);
saveas(fig,['cluster_prop.fig']);
close(fig);
% nodes table indicate matching in each cluster
nodtab=table(namesall',matchvec','VariableNames',{'nodes','match'});
writetable(nodtab,'nodetab.txt','Delimiter','\t');

% stack plot for clusters
runid=1;
loctab=est_tab(est_tab{:,'simulation'}==runid,:);
matchind=ppmmatch_ind_all(runid,:);
esttab_loc=loctab(matchind,:);
stackmat=[];
% unit convert
esttab_loc{:,1}=(esttab_loc{:,1}-para_add_list.conv_f(1))*para_add_list.conv_f(2);
groudtruth_onesamp_conv=groudtruth_onesamp;
groudtruth_onesamp_conv{:,1}=(groudtruth_onesamp_conv{:,1}-para_add_list.conv_f(1))*para_add_list.conv_f(2);
% the raw spectrum
tab_loc_arra=groudtruth_onesamp_conv{:,1:4};
sumsig=sin_mixture_simu(tab_loc_arra,timevec_sub_front,nan,'complex');
scalfactor=0.5;
sumsig(1)=sumsig(1)*scalfactor;
sumsig=[zeros([shifttimeadd,1]); sumsig];
spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,'temp');
stackmat=[stackmat; spec_new_sum{:,2}'];
%
estclusters=unique(graph_conn_group);
% the each cluster
ncluster=ncompound-1;
for groupi=1:ncluster%the first 4 clusters %(ncompound-1)%ignore DSS
  % the simulated peak groups truth
  groupind=find(groudtruth_onesamp_conv{:,'group'}==groupi);
  groundtruth_loc_arra=groudtruth_onesamp_conv{groupind,1:4};
  sumsig=sin_mixture_simu(groundtruth_loc_arra,timevec_sub_front,nan,'complex');
  scalfactor=0.5;
  sumsig(1)=sumsig(1)*scalfactor;
  sumsig=[zeros([shifttimeadd,1]); sumsig];
  spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,'temp');
  stackmat=[stackmat; spec_new_sum{:,2}'];
  % the matched cluster peaks
  groupind=find(graph_conn_group==estclusters(matchcluster(groupi)));
  esttab_loc_arra=esttab_loc{groupind,1:4};
  sumsig=sin_mixture_simu(esttab_loc_arra,timevec_sub_front,nan,'complex');
  scalfactor=0.5;
  sumsig(1)=sumsig(1)*scalfactor;
  sumsig=[zeros([shifttimeadd,1]); sumsig];
  spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig),imag(sumsig)),libdir,'temp');
  stackmat=[stackmat; spec_new_sum{:,2}'];
end
colorset=struct();
% colorset.rgb=flip([[0 0 0]; repmat([1 0 0; 0 0 1],[ncompound-1,1])],1);
% colorset.rgb=flip([[0 0 0]; repmat([1 0 0],[2,1]); repmat([0 1 0],[2,1]); repmat([0 0 1],[2,1]); repmat([1 1 0],[2,1])],1);
% colorset.categories=table(flip([{'sum'},{'cluster 1'},{'cluster 2'},{'cluster 3'},{'cluster 4'}]'));
% colorset.colorList=flip([0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0],1);
colorset.rgb=flip([[0 0 0]; repmat([1 0 0; 0 0 1],[ncluster,1])],1);
colorset.categories=table(flip([{'raw'},{'truth'},{'estimation'}]'));
colorset.colorList=flip([0 0 0; 1 0 0; 0 0 1],1);
% space
specvec=[0 1];
for clusteri=1:(ncluster-1)
  specvec=[specvec specvec(end)+9 specvec(end)+10];
end
specvec=[specvec specvec(end)+9]
%
region_loc=[6.8 8.6];
ppmvis_rang=sort(matchPPMs(region_loc,ppm));
ppmvis_ind=ppmvis_rang(1):ppmvis_rang(2);
stackmat=flip(stackmat,1);
% stackspec_time(stackmat(:,ppmvis_ind),ppm(ppmvis_ind)',0.0,300,['example_stack_plot'],'timeVect',specvec,'colors',colorset);
stackSpectra(stackmat(:,ppmvis_ind),ppm(ppmvis_ind)',0.0,300,['example_stack_plot'],'colors',colorset);
fig=gcf;
saveas(fig,['stack_clusters.fig']);
close all;
