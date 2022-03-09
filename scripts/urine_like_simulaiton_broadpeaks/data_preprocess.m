% construct the simulation dataset folder
%
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
% the path should be changed accordingly in the users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/urine_like_simulaiton_broadpeaks/'];
datadir=[paredir 'data/urine_fitting/'];
rundir=[projdir];
pipescriptdir=[localPaths.nmrdecomp_path 'scripts/urine_like_simualtion_phase/pipe_script/'];
preparapath=[datadir 'runid2_env_final.mat']; % decomposition of a urine dataset, to make the simulation more realistic. https://www.dropbox.com/s/pcbz1bgp9ss3kpb/runid2_env_final.mat?dl=0
preheadpath=[datadir 'test_trans.fid'];% a template fid file containing useful header information for NMRPipe. https://www.dropbox.com/s/1i0dixw4vasctwu/test_trans.fid?dl=0
% prenoisepath=[datadir 'temp.fid.txt'];%
cd(rundir);
mkdir('res')
cd([rundir '/res/']);
seedi=1;
rng(seedi);
% create the folder structure
%       ./res
%             ./deconv
%                   ./res (store deconvolution result from HPC)
%             ./nmrpipe_dir
%                   ./1...n (data to upload to HPC for decomposition)
%                   ./script (the nmrpipe scripts)
%             ./simu_data (store simulated spectra)
mkdir('./deconv');
mkdir('./deconv/res');
mkdir('./nmrpipe_dir');
mkdir('./nmrpipe_dir/script');
mkdir('./simu_data');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
nsample=20;
sampseq=[repmat(1,[10,1]); repmat(2,[10,1])];%10 without broad peaks + 10 with broad peaks
freq_res=600.132996;%MHZ
car_freq_ppm=4.763;%ppm
lambdarange=[1 15];
% load spectral region and add a broad peak
load(preparapath);
tabsumm_refine2=tabsumm_refine;%f, lambda, A, phi
remppm=[7.81 7.84; 7.55 7.555];
exist_ind=true(size(tabsumm_refine2,1),1);
for rowi=1:size(remppm,1);
  remf=(remppm(rowi,:)-para_add_list.conv_f(1))*para_add_list.conv_f(2);
  exist_ind=exist_ind & (tabsumm_refine2(:,1)<remf(1) | tabsumm_refine2(:,1)>remf(2));
end
tabsumm_refine2=tabsumm_refine2(exist_ind,:);
% refine signals compositons
tabsumm_refine2_sele=tabsumm_refine2;
[~,simu_ord]=sort(tabsumm_refine2_sele(:,1));
tabsumm_refine2_sele=tabsumm_refine2_sele(simu_ord,:);
nonneighb_ind=find(diff(tabsumm_refine2_sele(:,1))>=2);%remove peaks too close to each other
tabsumm_refine2_sele=tabsumm_refine2_sele(nonneighb_ind,:);
% randomize peak width lambda
tabsumm_refine2_sele(:,2)=rand([size(tabsumm_refine2_sele,1) 1])*(lambdarange(2)-lambdarange(1))+lambdarange(1);
% add DSS
tabsumm_refine2_sele=[tabsumm_refine2_sele; -car_freq_ppm*freq_res 4 0.1 0];
% add broad signals modifications
broadsignal=[1643,100,0.1,0; 1433,100,0.12,0; 1523,300,0.3,0];% ~7.5 ppm

% % plot check
% sumsig=sin_mixture_simu([tabsumm_refine2_sele; broadsignal],timevec_sub_front',1e-4,'complex');
% scalfactor=0.5;
% sumsig(1)=sumsig(1)*scalfactor;
% sumsig=[zeros([1,shifttimeadd]) sumsig];
% spec_new_sum=ft_pipe(table([1:length(sumsig)]',real(sumsig)',imag(sumsig)'),preheadpath,'temp');
% new_spec_vec2=spec_new_sum{:,2};
% fig=figure();
% plotr(ppm,new_spec_vec2,'LineWidth',2);%-mean(ft_ori_tab{:,2})
% hold on;
% stem(tabsumm_refine2_sele(:,1)/para_add_list.conv_f(2)+para_add_list.conv_f(1),tabsumm_refine2_sele(:,3));

ncompound=size(tabsumm_refine2_sele,1);
conc_range=[min(tabsumm_refine2_sele(:,3)) max(tabsumm_refine2_sele(:,3))];
sigma=1e-4;
% simulate spectra
fid_mat=[];
sepc_mat=[];
groundtruth_arra=[];
for sampi=1:nsample
  para_tab=tabsumm_refine2_sele;
  if sampseq(sampi)==2
    para_tab=[para_tab; broadsignal];
  end
  rng(rem(sampi,10));%so that there are pairs of samples with the same simulation except for the existence of broad peaks
  conc_vec=rand([ncompound-1 1])*(conc_range(2)-conc_range(1))+conc_range(1);
  conc_vec=[conc_vec; conc_range(2)];
  para_tab(1:length(conc_vec),3)=conc_vec;
  resvec=sin_mixture_simu(para_tab,timevec_sub_front',sigma,'complex');
  resvec(1)=resvec(1)*0.5;
  resvec=[zeros([1,shifttimeadd]) resvec];
  % spec_new_sum=ft_pipe(table([1:length(resvec)]',real(resvec)',imag(resvec)'),preheadpath,'temp');
  % sepc_mat=[sepc_mat; spec_new_sum{:,2}'];
  groundtruth_arra=[groundtruth_arra; [para_tab repmat(sampi,[size(para_tab,1),1])] ];
  fid_mat=[fid_mat; resvec];
end
groundtruth_tab=array2table(groundtruth_arra);
groundtruth_tab.Properties.VariableNames={'frequency','lambda','A','phase','simulation'};
timeind_vec=1:size(fid_mat,2);
% save as fid file
simusave=['simu_data/'];
for sampi=1:size(fid_mat,1)
  specvec=fid_mat(sampi,:);
  temptab=table(timeind_vec',real(specvec)',imag(specvec)');
  writetable(temptab,[simusave num2str(sampi) '.fid.txt'],'Delimiter',' ','WriteVariableNames',false);
end
% formulate the nmrpipe folder
nmrpipedir=['nmrpipe_dir/'];
for sampi=1:size(fid_mat,1)
  sampdir=num2str(sampi);
  mkdir([nmrpipedir sampdir]);
  mkdir([nmrpipedir sampdir '/script']);
  copyfile([simusave sampdir '.fid.txt'],[nmrpipedir sampdir '/temp.fid.txt']);
  copyfile([nmrpipedir 'script'],[nmrpipedir sampdir '/script']);
  copyfile(preheadpath,[nmrpipedir sampdir '/org.fid']);
  % copyfile(prenoisepath,[nmrpipedir sampdir '/noise.fid.txt']);
end

% nmrpipe based preprocess
nmrpipe_process('./nmrpipe_dir/',nsample,'prior');

% produce band by bucketing
ppm=[];
specmat=[];
for sampi=1:nsample
  samp=num2str(sampi);
  [value axes]=read_nmrp([nmrpipedir samp '/test.ft1']);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1';
  specmat=[specmat; value'];
end
%
% visulize the spctra
plotr(ppm,specmat);
specmat=flip(specmat,2);
ppm_r=flip(ppm);
%
[peaks.ints,peaks.shifts,peaks.params]=Peakpick1D(specmat,ppm_r,'max',0.6,'Complex');
hold off;

% use the selected paramters
size_bucket=0.002;
slackness=0.99;
[buckets]=optimize_optBucket(specmat,ppm_r,size_bucket,slackness);
clear('size_bucket','slackness')
% Filter out the bins with no peaks
[buckets]=filterBuckets_Peaks_opt(ppm_r,buckets,peaks);
%% click the figure
buckets=plotOptBucket_optResult(specmat,ppm_r,buckets,[7 8],[0 100]);
%%  Expand buckets
[buckets]=expandBucketBounds(buckets,specmat,ppm_r,'plotResult');
% get the bin ppm range
ppmrange=[-0.5 10];
binrange=buckets.optimized.expandedBuckets;
binrange=flip(binrange,2);
binrange=binrange(binrange(:,1)<ppmrange(2)&binrange(:,2)>ppmrange(1),:);%limited ppm searching range
% group every three bins
binrange=binrange';
binrange_vec=binrange(:);
%combine neighbor bins
groupind=[];
stepsize=6;
for i=1:floor(length(binrange_vec)/stepsize)
  num1=(i-1)*stepsize+1;
  num2=(i-1)*stepsize+stepsize;
  groupind=[groupind num1 num2];
end
binrange_comb_vec=binrange_vec(groupind);
binrange_comb=reshape(binrange_comb_vec,[2,length(binrange_comb_vec)/2])';%the range from high to low
% select targetted region
% ppmreg=[-0.5 10];
% binrange_comb=binrange_comb(max(binrange_comb,[],2)>ppmreg(1)&min(binrange_comb,[],2)<ppmreg(2),:);

% plot buckets
fig=figure(),
  hold on
  plotr(ppm_r,specmat);
  set(gca,'XDir','reverse')
  xlabel('Chemical Shift (ppm)')
  ylabel('Signal Intensity')
  set(gca, 'YTickLabel',[])
  title('Expanded Buckets')
  % Draw the new bin bounds
  highlightROIs(binrange_comb',max(specmat(:)),'color','r','edgeColor','k')
savefig(fig,'bin_select.fig');
close(fig);

writetable(table(binrange_comb),['binrange.txt'],'WriteVariableNames',false,'Delimiter','\t');
for sampi=1:nsample
  samp=num2str(sampi);
  copyfile(['./binrange.txt'],[nmrpipedir samp '/binrange.txt']);
end

binrange=binrange';
save('saved_simulation.mat');

nmrpipe_process('./nmrpipe_dir/',nsample,'band');

cd('./nmrpipe_dir');
!zip -r archive.zip *
cd('../');
