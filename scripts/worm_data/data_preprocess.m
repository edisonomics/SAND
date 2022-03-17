% construct the simulation dataset folder
% Different from other script, this script will be run on HPC as of the data size problem
% for worm data
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
paredir='/lustre2/scratch/yw44924/worm_nmr_fd/'
daradir=[paredir 'data/'];
rundir=[paredir];
pipescriptdir=[localPaths.nmrdecomp_path 'scripts/worm_nmr_fd/pipe_script/'];
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
mkdir('./worm_fid');
mkdir('./bin');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
% meta data table
meta_tab=readtable([datadir 'metadata.csv'],'Delimiter',',');
meta_tab=meta_tab(strcmp(meta_tab{:,'parameters'},'noesypr1d'),:);
deconvid=1:size(meta_tab,1);
meta_tab=addvars(meta_tab,deconvid','Before','set','NewVariableNames','deconvid');
parefold_pro='noesypr1d_batch';
folderslist={'noesypr1d_batch1and2' 'noesypr1d_batch3and4' 'noesypr1d_batch5and6'};
phasebase=[79.0 -8; 79.0 -8; 82.0 -8];%the basis phase for each folder
% adding phase correction information
meta_tab_new=[];
for folder_i=1:length(folderslist)
  folder=folderslist{folder_i};
  fileload=[datadir 'processed/' folder '/phase.csv'];
  phase_tab=readtable(fileload,'Delimiter',',');
  phase_tab.Properties.VariableNames={'file_order','p0','p1'};
  phase_tab{:,'p0'}=phase_tab{:,'p0'}+phasebase(folder_i,1);
  phase_tab{:,'p1'}=phase_tab{:,'p1'}+phasebase(folder_i,2);
  subfolder=repmat(folder,[size(phase_tab,1),1]);
  phase_tab=addvars(phase_tab,subfolder,'After','p1','NewVariableNames','subfolder');
  batch_info=cellfun(@str2num,strsplit(strrep(folder,'noesypr1d_batch',''),'and'));
  rowind=ismember(meta_tab{:,'batch'},batch_info);
  subtab=meta_tab(rowind,:);
  [~,orderind]=sort(subtab{:,'run_order'});
  file_order=1:size(subtab,1);
  subtab_ord=addvars(subtab(orderind,:),file_order','After','run_order','NewVariableNames','file_order');
  temptab=join(subtab_ord,phase_tab);
  meta_tab_new=[meta_tab_new; temptab];
end
meta_tab=meta_tab_new;
writetable(meta_tab,['metatab_all.txt'],'Delimiter','\t');

% format the fid folder
fidsave=['/res/worm_fid/'];
for rowi=1:size(meta_tab,1)
  locinfor=meta_tab(rowi,:);
  copyfile([datadir 'processed/' locinfor{:,'subfolder'} '/fid/' num2str(locinfor{:,'file_order'}) '.fid'],[fidsave num2str(locinfor{:,'deconvid'}) '.fid']);
end
% formulate the nmrpipe folder
nmrpipe=['nmrpipe_dir/'];
for rowi=1:size(meta_tab,1)
  locinfor=meta_tab(rowi,:);
  sampdir=num2str(locinfor{:,'deconvid'});
  mkdir([nmrpipe sampdir]);
  mkdir([nmrpipe sampdir '/script']);
  copyfile([fidsave sampdir '.fid'],[nmrpipe sampdir '/test.fid']);
  copyfile([nmrpipe 'script'],[nmrpipe sampdir '/script']);
  % modify for phase
  procshell=[nmrpipe sampdir '/script/proc.com'];
  % read lines
  % lines=readlines(procshell);
  lines={};
  f_id=fopen(procshell);
  while true
    tline=fgetl(f_id);
    if ~ischar(tline); break; end   %end of file
    % disp(tline);
    lines=[lines tline];
  end
  fclose(f_id);
  rowind=find(contains(lines,'-xELB'));
  for rowele=rowind
    lines{rowele}=['             ' '-xP0 ' num2str(locinfor{:,'p0'}) ' -xP1 ' num2str(locinfor{:,'p1'}) ' ' strtrim(lines{rowele})];
  end
  % output lines
  % cat(lines);
  f_id=fopen(procshell,'w');
  fprintf(f_id,'%s\n',lines{:});
  fclose(f_id);
end
% nmrpipe based preprocess
nmrpipe_process('./nmrpipe_dir/',nsample,'prior',shelladd);
% produce band by bucketing
ppm=[];
specmat=[];
for speci=1:size(meta_tab,1)
  sampdir=num2str(speci);
  [value axes]=read_nmrp([nmrpipe sampdir '/test.ft1']);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1';
  specmat=[specmat; value'];
end
% visulize the spctra
% plotr(ppm,specmat);
specmat=flip(specmat,2);
for speci=1:size(specmat,1)
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
