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
projdir=[paredir 'result_reproduce/krish_mixture/'];
datadir=[paredir 'data/'];
preresdirpath=[projdir 'res/deconv/res/res/'];
preheadpath=[projdir 'res/nmrpipe_dir/2/test_trans.fid'];
cd([projdir]);
%
load([projdir 'res/saved_simulation.mat']);
specppm=ppm_r;
%
nsample=33;
sampleseq=1:nsample;
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

% load decomposation estimation of different spectra
namelist={};
est_tab=[];%PPM, lambda, A, simulation_ind
ppmstr={};
ftstr={};
for sampi=sampleseq
  sample=sampleseq(sampi);
  samplestr=num2str(sample);
  load([preresdirpath samplestr '/runid' samplestr '_env_final.mat']);
  runtab=array2table(tabsumm_refine,'VariableNames',{'PPM','lambda','A','phase'});%f, lambda, A, phi
  nfeature=size(runtab,1);
  runtab{:,'PPM'}=runtab{:,'PPM'}/para_add_list.conv_f(2)+para_add_list.conv_f(1);
  simu=repmat(sampi,[nfeature,1]);
  runtab=[runtab table(simu)];
  namelist=[namelist; {repmat({'unknown'},[1,nfeature])}];
  est_tab=[est_tab; runtab];
  ppmstr{sampi}=ppm;
  ftstr{sampi}=ft_ori_tab;
end
est_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};

%intensity and integral based estimation
est_other_tab=[];
for sampi=sampleseq
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
  est_other_tab=[est_other_tab; runtab];
end
est_other_tab.Properties.VariableNames={'PPM','lambda','intensity','integral','phase','simulation'};
%
quan_str=struct();
quan_str.deconv=est_tab;
quan_str.intensity=est_other_tab(:,{'PPM','lambda','intensity','phase','simulation'});
quan_str.integral=est_other_tab(:,{'PPM','lambda','integral','phase','simulation'});
% select within interested ppm range
selerange=[0 13];% considered ppm range
remppms=[2.5 3.32];%the two highest peak related to DMSO and water
for type=fieldnames(quan_str)'
  type=type{1};
  tempdata=quan_str.(type);
  rowind_ppm=find(tempdata{:,'PPM'}>selerange(1)&tempdata{:,'PPM'}<selerange(2));
  tempdata=tempdata(rowind_ppm,:);
  remind=[];
  for remppm=remppms
    for sampi=1:nsample
      sampind=find(tempdata{:,'simulation'}==sampi);
      samp_tab=tempdata(sampind,:);
      regind=find(abs(samp_tab{:,'PPM'}-remppm)<0.05);
      [~,maxind]=max(samp_tab{regind,3});
      remind=[remind sampind(regind(maxind))];
    end
  end
  tempdata(remind,:)=[];
  % scale sum of quantification to 1
  for sampi=1:nsample
    sampind=find(tempdata{:,'simulation'}==sampi);
    tempdata{sampind,3}=tempdata{sampind,3}/sum(tempdata{sampind,3});
  end
  %
  quan_str.(type)=tempdata;
end
% calculate ground truth from the pure spectra
% use the peaks in pure spectra and relative concentration to calculate the expected peak intensity in each sample
groundtruth_ratio=readtable(['data/concentration_upd.txt']);
%
purespec_str=struct();
purespec_str.Ibuprofen=[1 12 23];
purespec_str.Prednesone=[11 22 33];
summ_str=struct();
groundtruth_str=struct();
for type=fieldnames(quan_str)'
  type=type{1};
  % ground truth peak list
  groundtruth_peak=[];
  for compd=fieldnames(purespec_str)'
    compd=compd{1};
    locseq=purespec_str.(compd);
    ppmlist_gt={};
    namelist_gt={};
    quantab=quan_str.(type);
    for locind=locseq
      locpeakind=quantab{:,'simulation'}==locind;
      locppm=quantab{locpeakind,'PPM'}';
      locA=quantab{locpeakind,3}';
      seleind=1:length(locppm);
      ppmlist_gt=[ppmlist_gt locppm(seleind)];
      namelist_gt=[namelist_gt {repmat({'unknown'},[1,length(seleind)])}];
    end
    [ppmmatch_ind_all ppmvec]=ppm_list_match(ppmlist_gt',namelist_gt,'^unknown$',deltapm_threshold);
    npeaks=length(ppmvec);
    A_arra=[];
    lambda_arra=[];
    for peaki=1:size(ppmmatch_ind_all,2)
      for locind_i=1:length(locseq)
        locind=locseq(locind_i);
        tempind=ppmmatch_ind_all(locind_i,peaki);
        loctab=quantab(quantab{:,'simulation'}==locind,:);
        A_arra(peaki,locind_i)=loctab{tempind,3};
        lambda_arra(peaki,locind_i)=loctab{tempind,'lambda'};
      end
    end
    compd_peak_tab=table(ppmvec',mean(lambda_arra,2),mean(A_arra,2),repmat(0,[npeaks,1]),repmat({compd},[npeaks,1]));
    seleind=1:npeaks;
    groundtruth_peak=[groundtruth_peak; compd_peak_tab(seleind,:)];
  end
  groundtruth_peak.Properties.VariableNames={'PPM','lambda','A_rela','phase','compounds'};
  % calculate ground truth A of peaks
  groundtruth_tab=[];
  conc_ref_compds={'Ibuprofen','Prednesone'};
  for sampi=sampleseq
    loctab=groundtruth_ratio(groundtruth_ratio{:,'ind'}==sampi,:);
    locconc=[loctab{:,'Ibuprofen'},loctab{:,'Prednesone'}];
    relconc=locconc./sum(locconc);
    relconc_exi=find(relconc>0);
    %
    for locind=relconc_exi
      matchcompd_mask=[];
      for rowi=1:size(groundtruth_peak,1)
        matchcompd_mask=[matchcompd_mask strcmp(groundtruth_peak{rowi,'compounds'},conc_ref_compds{locind})];
      end
      matchcompd_ind=find(matchcompd_mask);
      locpeak_tab=groundtruth_peak(matchcompd_ind,:);
      locpeak_tab{:,'A_rela'}=locpeak_tab{:,'A_rela'}*relconc(locind);
      locsavetabes=locpeak_tab(:,{'PPM','lambda','A_rela','phase','compounds'});
      locsavetabes=addvars(locsavetabes,repmat(sampi,[length(matchcompd_ind),1]),'After','phase');
      groundtruth_tab=[groundtruth_tab; locsavetabes];
    end
  end
  groundtruth_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation','compounds'};
  groundtruth_str.(type)=groundtruth_tab;
  % match estimation with ground truth
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
    subtab_true_match=subtab_true(ind_true,{'PPM','A','lambda','phase','simulation'});
    subtab_true_match.Properties.VariableNames={'PPM_true','A_true','lambda_true','phase_true','simulation'};
    summtab=[summtab; [subtab_est_match subtab_true_match]];
    rec_ratio=[rec_ratio size(subtab_est_match,1)/size(subtab_true,1)];
  end
  tempstr=struct();
  [~,sortind]=sort(summtab{:,'A_true'},'descend');
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
evalu_tab=cell2table(cell(0,5),'VariableNames',{'rel_mse','mse','corxy','k', 'quan_method'});
for methele=fieldnames(evalu_str)'
  methele=methele{1};
  loctab=struct2table(evalu_str.(methele));
  loctab=[loctab(:,{'rel_mse','mse','corxy','k'}) table({methele},'VariableNames',{'quan_method'})];
  evalu_tab=[evalu_tab; loctab];
end
save('evaluation.mat','evalu_str','evalu_tab');
writetable(evalu_tab,'stat_tab.txt');

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
est_tab_deconv=quan_str.deconv;
est_tab_dec_subset={};
groundtruth_tab_deconv=groundtruth_str.deconv;
samp_cor_seq=setdiff(1:nsample,[purespec_str.Ibuprofen purespec_str.Prednesone]);%only samples that contains the two compounds
for runid=samp_cor_seq
  loctab=est_tab_deconv(est_tab_deconv{:,'simulation'}==runid,:);
  ppmlist=[ppmlist loctab{:,'PPM'}'];
  namelist=[namelist {repmat({'unknown'},[1,size(loctab,1)])}];
  intlist=[intlist loctab{:,'A'}];
  est_tab_dec_subset=[est_tab_dec_subset {loctab}];
end
[ppmmatch_ind_all ppmvec]=ppm_list_match(ppmlist',namelist,'^unknown$',deltapm_threshold);
mat_reshape_all=[];
for isample=1:size(ppmmatch_ind_all,1)
  intenvec=intlist{isample};
  ppmind=ppmmatch_ind_all(isample,:);
  intenarray=intenvec(ppmind);
  mat_reshape_all=[mat_reshape_all; intenarray'];
  est_tab_dec_subset{isample}=est_tab_dec_subset{isample}(ppmind,:);
end
%% match with ground truth peak sets
subtab_true=groundtruth_tab_deconv(groundtruth_tab_deconv{:,'simulation'}==samp_cor_seq(1),:);
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
for runid=1:length(est_tab_dec_subset)
  est_tab_dec_subset{runid}=est_tab_dec_subset{runid}(ind_est,:);
end
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
cormat_thres=cor_intensity>=0.5;%thres_display_vec_inten;
G=graph(cormat_thres,'omitselfloops');
% plot(G);
[graph_conn_group]=conncomp(G);
groudtruth_onesamp=groundtruth_tab_deconv(groundtruth_tab_deconv{:,'simulation'}==samp_cor_seq(1),:);
est_tab_dec_onesamp=est_tab_dec_subset{1};
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
matchratio=[];
visregion=[0 2.5; 4 8];
for compd=fieldnames(purespec_str)'
  compd=compd{1};
  for regioni=1:size(visregion,1)
    locregion=visregion(regioni,:);
    compdind=find(strcmp(groudtruth_onesamp{:,'compounds'},compd));
    groupsize=length(compdind);
    matchedind=find(ismember(ind_true,compdind));
    coll_clust=[];
    g_group_unique=unique(graph_conn_group);
    for estclut=g_group_unique
      groupind=find(graph_conn_group==estclut);
      memmatch=find(ismember(ind_est(matchedind),groupind));
      coll_clust=[coll_clust length(memmatch)];
    end
    [maxmatch,maxind]=max(coll_clust);
    matchratio=[matchratio maxmatch/groupsize];
    % recovered spectra for each clusters
    clusind=intersect(ind_true,compdind);
    tabsumm_refine2_sele=groudtruth_onesamp{clusind,1:4};
    tabsumm_refine2_sele(:,1)=(tabsumm_refine2_sele(:,1)-para_add_list.conv_f(1))*para_add_list.conv_f(2);
    sumsig=sin_mixture_simu(tabsumm_refine2_sele,timevec_sub_front',0.0,'complex');
    scalfactor=0.5;
    sumsig(1)=sumsig(1)*scalfactor;
    sumsig=[zeros([1,shifttimeadd]) sumsig];
    spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig)',imag(sumsig)'),preheadpath,'temp');
    new_spec_vec2=spec_new_sum{:,2};
    % pure sample spectra
    sampplot_inds=purespec_str.(compd);
    ft_ori_mat=[];
    for puresampind=sampplot_inds
      ppmori=ppmstr{puresampind};
      ft_ori_mat=[ft_ori_mat ftstr{puresampind}{:,2}];
    end
    newinten=ft_ori_mat(:,1);
    ppmignore=[2.2 3.5];
    igreg=sort(matchPPMs(ppmignore,ppmori));
    newinten_reg=newinten([1:igreg(1) igreg(2):length(newinten)]);
    int_conv_fac=max(new_spec_vec2)/max(newinten_reg);
    % local ppm region
    reg1=sort(matchPPMs(locregion,ppmori));
    seq1=reg1(1):reg1(2);
    reg2=sort(matchPPMs(locregion,ppm));
    seq2=reg2(1):reg2(2);
    %
    fig=figure();
    plotr(ppmori(seq1),ft_ori_mat(seq1,:)*int_conv_fac,'LineWidth',2,'Color',[0.8 0.8 0.8]);
    hold on;
    plotr(ppm(seq2),new_spec_vec2(seq2),'LineWidth',2,'Color','r');%-mean(ft_ori_tab{:,2})
    % legend('raw','estimation');
    title(compd);
    savefig(fig,[compd '_' num2str(regioni) '_cluster.fig']);
    close(fig);
  end
end
