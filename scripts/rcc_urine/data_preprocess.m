% construct the dataset folder
% Different from other script, this script will be run on HPC as of the data size problem
% for urine data
close all;
clear all;
shelladd='source /PATHTO/.cshrc';
%% Set your toolbox paths; functions imported from these directories:
% rmpath(genpath('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/matalb.lib/plstoolbox/PLS_Toolbox_881/'));
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/PATHTO/Edison_Lab_Shared_Metabolomics_UGA/';
% NMR decompositon program found @ https://github.com/edisonomics/SAND
localPaths.nmrdecomp_path='/PATHTO/SAND/';
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.nmrdecomp_path));
pause(1),clc
% the path should be changed accordingly in the users' computer
paredir='/PATHTO/rcc_urine/'
datadir=[paredir 'data/'];
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
mkdir('./deconv/temp');
mkdir('./deconv/data');
mkdir('./nmrpipe_dir');
mkdir('./nmrpipe_dir/script');
mkdir('./urine_fid');
mkdir('./bin');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
%
meta_tab=readtable([datadir 'meta_data_1/sample_info_all_groups.xlsx']);
nsample=size(meta_tab,1);
rawdatadir=[datadir 'BR1_BIF_5_noesypr1d/'];
% get the local folder list
cd(rawdatadir);
c=strread(ls,'%s');
c=sort(cellfun(@str2num,c(find(~(cellfun(@isempty,regexp(c,'^\d+$')))))));
datafds=num2str(c);
cd([rundir '/res/']);
%
meta_tab=addvars(meta_tab,datafds,'After','Yvec');
% formulate the nmrpipe folder
nmrpipe=['nmrpipe_dir/'];
for rowi=1:nsample
  sampdir=num2str(rowi);
  oridir=strtrim(meta_tab{rowi,'datafds'});
  mkdir([nmrpipe sampdir]);
  mkdir([nmrpipe sampdir '/script']);
  copyfile([rawdatadir oridir '/*'],[nmrpipe sampdir '/']);
  copyfile([nmrpipe 'script'],[nmrpipe sampdir '/script']);
end
% nmrpipe based preprocess
nmrpipe_process('./nmrpipe_dir/',nsample,'prior',shelladd);

% produce band by bucketing
ppm=[];
specmat=[];
for speci=1:nsample
  sampdir=num2str(speci);
  [value axes]=read_nmrp([nmrpipe sampdir '/test.ft1']);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1';
  specmat=[specmat; value'];
end
% visulize the spctra
% plotr(ppm,specmat);
specmat=flip(specmat,2);
for speci=1:nsample
  tempspec_mat=specmat(speci,:);
  tempspec_mat=[tempspec_mat; tempspec_mat];%because the peak picking cannot work with one spectra
  ppm_r=flip(ppm);
  [peaks.ints,peaks.shifts,peaks.params]=Peakpick1D(tempspec_mat,ppm_r,'max',0.6,'Complex');
  hold off;
  % % parameter searching
  % size_bucket=0.002:0.002:0.01;
  % slackness=0.45:0.13:0.99;
  % [buckets]=optimize_optBucket(tempspec_mat,ppm_r,size_bucket,slackness);
  % clear('size_bucket','slackness')
  % [buckets]=filterBuckets_Peaks_opt(ppm_r,buckets,peaks);
  % buckets=plotOptBucket_optResult(tempspec_mat,ppm_r,buckets,[3.7 3.9],[0 150]);

  % use the selected paramters
  size_bucket=0.002;
  slackness=0.99;
  [buckets]=optimize_optBucket(tempspec_mat,ppm_r,size_bucket,slackness);
  clear('size_bucket','slackness')
  % Filter out the bins with no peaks
  [buckets]=filterBuckets_Peaks_opt(ppm_r,buckets,peaks);
  close all;
  %% automatically click the figure
  buckets=plotOptBucket_optResult_auto(tempspec_mat,ppm_r,buckets,[3.5 4.0],[0 150]);
  close all;
  %%  Expand buckets
  [buckets]=expandBucketBounds(buckets,tempspec_mat,ppm_r,'plotResult');
  % [buckets2]=expandBucketBounds(buckets,tempspec_mat([3 44 49 92 97 144],:),ppm_r,'plotResult');
  % [buckets2]=expandBucketBounds(buckets,tempspec_mat([1:10],:),ppm_r,'plotResult');
  % get the bin ppm range
  ppmrange=[-0.5 10];
  binrange=buckets.optimized.expandedBuckets;
  binrange=flip(binrange,2);
  binrange=binrange(binrange(:,1)<ppmrange(2)&binrange(:,2)>ppmrange(1),:);%limited ppm searching range
  % group every three bins
  binrange=binrange';
  binrange_vec=binrange(:);
  %combine each neighbor 3 range
  groupind=[];
  stepsize=6;
  for i=1:floor(length(binrange_vec)/stepsize)
    num1=(i-1)*stepsize+1;
    num2=(i-1)*stepsize+stepsize;
    groupind=[groupind num1 num2];
  end
  binrange_comb_vec=binrange_vec(groupind);
  binrange_comb=reshape(binrange_comb_vec,[2,length(binrange_comb_vec)/2])';%the range from high to low
  % plot buckets
  fig=figure(),
    hold on
    plotr(ppm_r,tempspec_mat);
    set(gca,'XDir','reverse')
    xlabel('Chemical Shift (ppm)')
    ylabel('Signal Intensity')
    set(gca, 'YTickLabel',[])
    title('Expanded Buckets - Lowest points method')
    % Draw the new bin bounds
    highlightROIs(binrange_comb',max(tempspec_mat(:)),'color','r','edgeColor','k')
  saveas(fig,['./bin/' num2str(speci) '_bin.fig']);
  
  writetable(table(binrange_comb),['./bin/' num2str(speci) '_binrange.txt'],'WriteVariableNames',false,'Delimiter','\t');
  sampdir=num2str(speci);
  copyfile(['./bin/' sampdir '_binrange.txt'],[nmrpipe sampdir '/binrange.txt']);
  close all hidden;
end
binrange=binrange';
save('saved_preprocessing.mat');

nmrpipe_process('./nmrpipe_dir/',nsample,'band',shelladd);

copyfile(['../runtab.txt'],'./deconv/');

!mv ./nmrpipe_dir/* ./deconv/data/
