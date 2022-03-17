% construct the simulation dataset folder
% Different from other script, this script will be run on HPC as of the data size problem
% for bif's urine data
% 1 = Control
% 2 = Pure Clear Cell RCC
% 3 = Cancer (other) - unclassified
% 4 = Papillary RCC
% 5 = Chromophobe Carcinoma
% 6 = Sarcomatoid
% 7 = Clear Cell Papillary
% 8 = External Control
% 9 = Internal Control
close all;
clear all;
shelladd='source /home/yw44924/.cshrc';
%% Set your toolbox paths; functions imported from these directories:
% rmpath(genpath('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/matalb.lib/plstoolbox/PLS_Toolbox_881/'));
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/home/yw44924/metabolic_toolbox_open/Edison_Lab_Shared_Metabolomics_UGA/';
% NMR decompositon program found @ https://github.com/mikeaalv/NMR_time_domain_decomposition
localPaths.nmrdecomp_path='/home/yw44924/nmr_spec_decomp/NMR_time_domain_decomposition/';
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.nmrdecomp_path));
pause(1),clc
% the path should be changed accordingly in the users' computer
paredir='/lustre2/scratch/yw44924/bif_urine/'
daradir=[paredir 'data/'];
rundir=[paredir];
pipescriptdir=[localPaths.nmrdecomp_path 'scripts/rcc_urine/pipe_script/'];
%
cd(rundir);
mkdir('res')
cd([rundir '/res/']);
seedi=1;
rng(seedi);
% create the folder structure
%       ./data
%       ./res
%             ./deconv
%                   ./res (the file to be downloaded)
%             ./nmrpipe_dir
%                   ./1...n (data to upload to HPC for decomposition)
%                   ./script (the nmrpipe scripts)
%             ./ft_data (store simulated spectra)
%             a few files to be downloaded
mkdir('./deconv');
mkdir('./deconv/res');
mkdir('./nmrpipe_dir');
mkdir('./nmrpipe_dir/script');
mkdir('./ft_data');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
% load and setup spectra
spectra=Load1D([daradir 'BR1_BIF_5_noesypr1d'],'bruker');
[X,ppm]=Setup1D(spectra);
metatab=readtable([daradir 'meta_data_1/sample_info_all_groups.xlsx']);
Yvec_meta=metatab.Yvec;
% reference
%% Referencing to DSS (manual)
%spectra = ref_spectra_manual(spectra, [-0.02, 0.01], 0)
%% Reference spectra
spectra1=ref_spectra_max(spectra,-0.02)
% whichLine()
[X,ppm]=Setup1D(spectra1);
displaypeak1D(X,ppm,0,Yvec_meta);
fig=gcf;
saveas(fig,['spectra_referenced.fig'])
close all;
%% Removing ends and water regions
% Spectra ends are removed together with water region (4.683-4.893)
[~,k1]=min(abs(ppm--0.5));
XR=X(:,k1:end);
ppmR=ppm(k1:end);
[~,k2]=min(abs(ppmR-10.0));
XR1=XR(:,1:k2);
ppmR=ppmR(1:k2);
XR2=remove_region(XR1,ppmR,4.893,4.683);
% Alignment using CCOW
% XAL=XR2;
XAL=guide_align1D(XR2,ppmR,'spearman','CCOW');
displaypeak1D(XAL,ppmR,0,Yvec_meta);
fig=gcf;
saveas(fig,['spectra_aligned.fig'])
close all;
ppm_Xaligned=vertcat(ppmR,XAL);
csvwrite('class_aligned.csv',ppm_Xaligned')
% Normalization using probabilistic quotient normalization (PQN) method
XALN=normalize(XAL,ppmR,'PQN');
normcheck(XALN)
fig=gcf;
saveas(fig,['normalcheck_after.fig'])
close all;
displaypeak1D(XALN,ppmR,0,Yvec_meta);
fig=gcf;
saveas(fig,['spectra_normalized.fig'])
close all;
ppm_XALN_PQN=vertcat(ppmR,XALN);
csvwrite('class_PQN.csv',ppm_XALN_PQN')
% Remove a Spectra.
%'SSO54' Papillary RCC sample was removed from furthur analysis because of poor spectral quality, and becaue it prevented the quantification of several metabolites.
% This sample is line 269 - pre-removal.
XALN=XALN([1:268, 270:end],:);
%% Delete SSO54 (269) from T table.
metatab=readtable([daradir 'meta_data_1/sample_removedSS054.xlsx']);
Yvec_meta=metatab.Yvec;
%
nsample=size(XALN,1);
% save the aligned spectra
% displaypeak1D(XALN,ppmR,0,Yvec);
% plot(ppmR,XALN)
presave=['ft_data/'];
for sampi=1:size(XALN,1)
  specvec=XALN(sampi,:);
  temptab=table(flip(ppmR)',flip(specvec)');%as nmrpipe assume one direction of reading spectra
  writetable(temptab,[presave num2str(sampi) '.ft.txt'],'Delimiter',' ','WriteVariableNames',false);
end
% formulate the nmrpipe folder
nmrpipedir=['nmrpipe_dir/'];
for sampi=1:size(XALN,1)
  sampdir=num2str(sampi);
  mkdir([nmrpipedir sampdir]);
  mkdir([nmrpipedir sampdir '/script']);
  copyfile([presave sampdir '.ft.txt'],[nmrpipedir sampdir '/temp.ft.txt']);
  copyfile([nmrpipedir 'script'],[nmrpipedir sampdir '/script']);
  % copyfile(preheadpath,[nmrpipedir sampdir '/org.fid']);
  % copyfile(prenoisepath,[nmrpipedir sampdir '/noise.fid.txt']);
end
% PCA plot
%% Scale using 'logoff'
XALNS2=scale(XALN,'logoff');
varcheck(XALNS2)
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
% produce band by bucketing
[peaks.ints,peaks.shifts,peaks.params]=Peakpick1D(XALN,ppmR,'max',0.6,'Complex');
hold off;
% use the selected paramters
size_bucket=0.002;
slackness=0.99;
[buckets]=optimize_optBucket(XALN,ppmR,size_bucket,slackness);
clear('size_bucket','slackness')
% Filter out the bins with no peaks
[buckets]=filterBuckets_Peaks_opt(ppmR,buckets,peaks);
%% click the figure
close all;
buckets=plotOptBucket_optResult_auto(XALN,ppmR,buckets,[7 8],[0 0.1]);
%%  Expand buckets
[buckets]=expandBucketBounds(buckets,XALN,ppmR,'plotResult');
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
  plotr(ppmR,XALN);
  set(gca,'XDir','reverse')
  xlabel('Chemical Shift (ppm)')
  ylabel('Signal Intensity')
  set(gca, 'YTickLabel',[])
  title('Expanded Buckets')
  % Draw the new bin bounds
  highlightROIs(binrange_comb',max(XALN(:)),'color','r','edgeColor','k')
savefig(fig,'bin_select.fig');
close(fig);

writetable(table(binrange_comb),['binrange.txt'],'WriteVariableNames',false,'Delimiter','\t');
for sampi=1:nsample
  samp=num2str(sampi);
  copyfile(['./binrange.txt'],[nmrpipedir samp '/binrange.txt']);
end

binrange=binrange';
save('saved_preprocessing.mat');

nmrpipe_process('./nmrpipe_dir/',nsample,'band',shelladd);

% cd('./nmrpipe_dir');
% !zip -r archive.zip *
% cd('../');

% create the running table
runtab=readtable([daradir 'runtab_template.txt']);
runidseq=cellfun(@(x) {num2str(x)},num2cell([1:nsample]))';
lambdamax=repmat(unique(runtab{:,'lambdamax'}),[nsample,1]);
newprop=repmat(unique(runtab{:,'newprop'}),[nsample,1]);
thres_digit=repmat(unique(runtab{:,'thres_digit'}),[nsample,1]);
multi_replicate=repmat(unique(runtab{:,'multi_replicate'}),[nsample,1]);
dataset=cellfun(@(x) {['../nmrpipe_dir/' x '/']},runidseq);
runtab_new=table(runidseq,lambdamax,newprop,thres_digit,multi_replicate,dataset);
writetable(runtab_new,['runtab.txt'],'Delimiter','\t','WriteVariableNames',false);
copyfile(['runtab.txt'],'./deconv/');
