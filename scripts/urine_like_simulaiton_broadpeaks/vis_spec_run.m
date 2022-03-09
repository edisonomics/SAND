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
projdir=[paredir 'result_reproduce/urine_like_simulaiton_broadpeaks/'];
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
groundtruth_tab.Properties.VariableNames={'PPM','lambda','A','phase','simulation'};

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
    subtab_true_match=subtab_true(ind_true,{'PPM','A','lambda','phase','simulation'});
    subtab_true_match.Properties.VariableNames={'PPM_true','A_true','lambda_true','phase_true','simulation'};
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
save('evaluation.mat','evalu_str');
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
  title([type ' lambda ']);
saveas(h,['scatter_simulation.' type '_lambda.fig']);
close(h);
corr(summtab{:,'lambda_true'},summtab{:,'lambda_est'})
% phi estimation
unique(summ_str.deconv.summtab{:,'phase_est'})
