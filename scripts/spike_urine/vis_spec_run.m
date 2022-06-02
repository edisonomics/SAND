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
preresdirpath=[projdir 'res/res/'];
cd([projdir]);
load([projdir 'res/saved_preprocessing.mat']);
specppm=ppm_r;
%
sampleseq=1:nsample;
ppmrange_dss=[-0.1 0.1];
deltapm_threshold=0.01;%distance threshold for peak matching
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

stackSpectra(ft_mat,ppm,0,50,'all stack')
fig=gcf;
saveas(fig,['stack_all_sample.fig']);
close all;
% load decomposation estimation of different spectra
namelist={};
est_tab=[];%PPM, lambda, A, simulation_ind
ppmstr={};
ftstr={};
for sampi=sampleseq
  samplestr=num2str(sampi);
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
  samplestr=num2str(sampi);
  load([preresdirpath samplestr '/runid' samplestr '_env_final.mat']);
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
%
selerange=[-0.5 9];% considered ppm range
dssrange=[-0.5 0.5];
% ratios
spike_ratio=[0,1,2,3,4];%ratio for spike-in compounds
dss_ratio=(1/9*0.06+0.4*spike_ratio)/(1/9*0.06+0.4);%ratio for dss
% convert factors
gt_ratio=spike_ratio.*spike_ratio./dss_ratio;
est_ratio=spike_ratio;
% group ind
mixseq=[2:5];%mixutres
mixseq=[4:5];%remove the low concentration sample
purespec_str.ref=[6];%reference sample
remreg=[4.66 4.95];% the water region to remove
for type=fieldnames(quan_str)'
  type=type{1};
  tempdata=quan_str.(type);
  rowind_ppm=find(tempdata{:,'PPM'}>selerange(1)&tempdata{:,'PPM'}<selerange(2));
  tempdata=tempdata(rowind_ppm,:);
  %
  remind=[];
  for sampi=1:nsample
    sampind=find(tempdata{:,'simulation'}==sampi);
    samp_tab=tempdata(sampind,:);
    regind=find(samp_tab{:,'PPM'}>remreg(1)&samp_tab{:,'PPM'}<remreg(2));
    remind=[remind sampind(regind)'];
  end
  tempdata(remind,:)=[];
  % scale by dss peak and apply convert ratios
  dssindvec=[];
  for sampi=1:nsample
    sampind=find(tempdata{:,'simulation'}==sampi);
    loctab=tempdata(sampind,:);
    dssrangind=find(loctab{:,1}>dssrange(1)&loctab{:,1}<dssrange(2));
    [a_dss, dssind]=max(loctab{dssrangind,3});
    loctab{:,3}=loctab{:,3}/a_dss;
    if ismember(sampi,mixseq)
      loctab{:,3}=loctab{:,3}*est_ratio(sampi);
    end
    dssindvec=[dssindvec sampind(dssrangind(dssind))];
    tempdata(sampind,:)=loctab;
  end
  %
  % tempdata(dssindvec,:)=[];
  quan_str.(type)=tempdata;
end
% calculate ground truth from the reference spectra
% match ppms
reg_ref=[1.4705 1.4826 3.0521 3.0644 3.2223 3.236 3.2379 3.2511 3.3768 3.3876 3.3925 3.4027 3.4081 3.4189 3.4423 3.4462 3.4521 3.4555 3.4663 3.4819 3.4971 3.5147 3.521 3.5313 3.5377 3.5963 3.6037 3.8736 3.877 3.8941 3.898 3.6305 4.6437 5.2216 5.228 6.8948 6.9095 7.3491 7.3447 7.8586 7.863];
reg_spike=[1.4743 1.486 3.0536 3.0663 3.2296 3.2433 3.2453 3.2585 3.3856 3.3964 3.4013 3.412 3.4174 3.4281 3.4477 3.4516 3.4575 3.4609 3.4761 3.4912 3.5064 3.5245 3.5308 3.5411 3.547 3.6066 3.6139 3.8824 3.8858 3.9029 3.9068 4.6398 3.653 5.2309 5.2373 6.933 6.9476 7.3784 7.374 7.8737 7.8781];
%
summ_str=struct();
groundtruth_str=struct();
ground_peak_n=200;
for type=fieldnames(quan_str)'
  type=type{1};
  % ground truth peak list
  groundtruth_peak=[];
  locseq=purespec_str.ref;
  quantab=quan_str.(type);
  locpeakind=quantab{:,'simulation'}==locseq;
  groundtruth_peak=[groundtruth_peak; [quantab(locpeakind,1:4) table(repmat({'unknown'},[length(find(locpeakind)),1]))]];
  groundtruth_peak.Properties.VariableNames={'PPM','lambda','A_rela','phase','compounds'};
  [~,indsort]=sort(groundtruth_peak{:,3},'descend');
  groundtruth_peak=groundtruth_peak(indsort(1:ground_peak_n),:);
  % calculate ground truth A of peaks in each sample
  groundtruth_tab=[];
  for sampi=mixseq
    locpeak_tab=groundtruth_peak;
    locpeak_tab{:,'A_rela'}=locpeak_tab{:,'A_rela'}*gt_ratio(sampi);
    locsavetabes=locpeak_tab(:,{'PPM','lambda','A_rela','phase','compounds'});
    locsavetabes=addvars(locsavetabes,repmat(sampi,[size(locsavetabes,1),1]),'After','phase');
    groundtruth_tab=[groundtruth_tab; locsavetabes];
  end
  groundtruth_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation','compounds'};
  groundtruth_str.(type)=groundtruth_tab;
  % match estimation with ground truth
  summtab=[];
  rec_ratio=[];
  est_tab_temp=quan_str.(type);
  est_tab_temp.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};
  for simui=mixseq
    subtab_est=est_tab_temp(est_tab_temp{:,'simulation'}==simui,:);
    subtab_true=groundtruth_tab(groundtruth_tab{:,'simulation'}==simui,:);
    % two step matching
    % peaks that moves
    ind_true_match=[];
    ind_est_match=[];
    for regi=1:size(reg_ref,2)
      [~,truesetind]=min(abs(subtab_true{:,'PPM'}-reg_ref(regi)));
      [~,estsetind]=min(abs(subtab_est{:,'PPM'}-reg_spike(regi)));
      ind_est_match=[ind_est_match; estsetind];
      ind_true_match=[ind_true_match; truesetind];
    end
    % peaks that doesn't move
    subtab_est_nomove=subtab_est;
    subtab_true_nomove=subtab_true;
    subtab_est_nomove(ind_est_match,:)=[];
    subtab_true_nomove(ind_true_match,:)=[];
    distmat=abs(pdist2(subtab_est_nomove{:,'PPM'},subtab_true_nomove{:,'PPM'}));
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
    subtab_est_match=[subtab_est(ind_est_match,{'PPM','A','lambda','phase'}); subtab_est_nomove(ind_est,{'PPM','A','lambda','phase'})];
    subtab_est_match.Properties.VariableNames={'PPM_est','A_est','lambda_est','phase_est'};
    subtab_true_match=[subtab_true(ind_true_match,{'PPM','A','lambda','phase','simulation'}); subtab_true_nomove(ind_true,{'PPM','A','lambda','phase','simulation'})];
    subtab_true_match.Properties.VariableNames={'PPM_true','A_true','lambda_true','phase_true','simulation'};
    summtab=[summtab; [subtab_est_match subtab_true_match]];
    rec_ratio=[rec_ratio size(subtab_est_match,1)/size(subtab_true,1)];
  end
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
  for simui=mixseq
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
    gscatter(log10(summtab{:,'A_true'}),log10(summtab{:,'A_est'}),summtab{:,'simulation'},[],[],[20]);
    xlabel('ground truth');
    ylabel('estimation');
    title([type ' correlation ' num2str(evalu.corxy), ' mse ' num2str(evalu.mse) ' k ' num2str(evalu.k)]);
  saveas(h,['scatter_simulation.' type 'log10.fig']);
  close(h);
end

% example regions


% correlation network
%% match peaks from different samples
ppmlist={};
namelist={};
intlist={};
est_tab_deconv=quan_str.deconv;
est_tab_dec_subset={};
groundtruth_tab_deconv=groundtruth_str.deconv;
for runid=mixseq
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
subtab_true=groundtruth_tab_deconv(groundtruth_tab_deconv{:,'simulation'}==mixseq(1),:);
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
