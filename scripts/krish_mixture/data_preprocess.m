% construct the simulation dataset folder
%
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
% the path should be changed accordingly in the users' computer
paredir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/publicaiton_spec_decomp/'
projdir=[paredir 'result_reproduce/krish_mixture/'];
datadir=[projdir 'data/'];
rundir=[projdir];
pipescriptdir=[localPaths.nmrdecomp_path 'scripts/krish_mixture/pipe_script/'];
%
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
%             ./raw_fid (store simulated spectra)
mkdir('./deconv');
mkdir('./deconv/res');
mkdir('./nmrpipe_dir');
mkdir('./nmrpipe_dir/script');
mkdir('./raw_fid');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
% format the fid folder
nsample=33;
fidsave=[paredir '/res/raw_fid/'];
for rowi=1:nsample
  sampdir=num2str(rowi);
  mkdir([fidsave sampdir]);
  copyfile([datadir 'fid_raw_renamed/' sampdir '/*'],[fidsave sampdir '/']);
end
% formulate the nmrpipe folder
nmrpipe=['nmrpipe_dir/'];
for rowi=1:nsample
  sampdir=num2str(rowi);
  mkdir([nmrpipe sampdir]);
  mkdir([nmrpipe sampdir '/script']);
  copyfile([fidsave sampdir '/*'],[nmrpipe sampdir '/']);
  copyfile([nmrpipe 'script'],[nmrpipe sampdir '/script']);
end
% nmrpipe based preprocess
nmrpipe_process('./nmrpipe_dir/',nsample,'prior');

% produce band by bucketing
nmrpipedir='nmrpipe_dir/';
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
buckets=plotOptBucket_optResult(specmat,ppm_r,buckets,[7 8],[0 200]);
%%  Expand buckets
[buckets]=expandBucketBounds(buckets,specmat,ppm_r,'plotResult');
% get the bin ppm range
ppmrange=[-0.5 13];
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
