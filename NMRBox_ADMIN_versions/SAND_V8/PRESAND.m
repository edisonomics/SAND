%YW wrote the base
%preprocess SAND V7: OS modified some 08-02-2022
%This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.

%this section is only for debugging: make sure to run it under the main
%root
% setenv('SANDPATH','/home/nmrbox/osanati/read-only/test1')
% setenv('ORIGIN','/home/nmrbox/osanati/read-only/usert1')
% setenv('USERF',[getenv('ORIGIN'), '/try1'])
% setenv('USERL',[getenv('ORIGIN'), '/SAND_last_run'])
% setenv('TDD',[getenv('SANDPATH'),'/src']) 
% setenv('ELSM',[getenv('SANDPATH'),'/Edison_Lab_Shared_Metabolomics_UGA'])
% setenv('NCPUS','5')
% setenv('WDIR', getenv('USERL'))
% setenv('TEMPDIR',[getenv('USERL'),'temp'])
% setenv('LOGDIR',[getenv('USERL'),'/log'])
% setenv('NSAMPLE','1')
% setenv('NSIG','7')  %max number of peaks in each bin'
% setenv('NTRY','2000') %number of iteration')
% setenv('TEMP','80') %temparature in MCMC
% setenv('PREP','1')
% setenv('LPPM','1.1')
% setenv('UPPM','1.16')
% setenv('GBIN','3')

close all;
clear all;
shelladd=':';     
localPaths.public_toolbox=getenv('ELSM')  %Edison_Lab_Shared_Metabolomics_UGA'
localPaths.nmrdecomp_path=getenv('TDD') 
addpath(genpath(getenv('ELSM')));
addpath(genpath(getenv('TDD')));
pause(1),clc

lastdir=getenv('USERL'); %where to write last_results rundir: USERL=$ORIGIN/SAND_last_run
lastdir=[lastdir,'/']
libdir=getenv('SANDPATH'); %where scripts are
datadir=getenv('USERF'); %where the user NMR folders data is. datadir=[paredir '/data/'];
datadir=[datadir '/'];

datadir2=datadir; %should be changed for unprocessed data

pipescriptdir=[libdir '/pipe_scripts/'];
%
cd(lastdir);
system('rm -f -r ./res');
mkdir('res')
cd([lastdir 'res/']);
seedi=1;
rng(seedi);
% create the folder structure
%       ./data
%       ./res
%             ./nmrpipe_dir
%                   ./1...n (data to upload to HPC for decomposition)
%                   ./script (the nmrpipe scripts)
mkdir('./nmrpipe_dir');
mkdir('./nmrpipe_dir/script');
mkdir('./bin');
copyfile([pipescriptdir '*'],'./nmrpipe_dir/script/');
%
 
% get the local folder list and copy all data to last run folder with adjusted number: assume only folders are there
if (str2num(getenv('PREP')) == 1) 
	cd(datadir);
	rawfolders=strsplit(strtrim(ls));
	nsample=size(rawfolders,2);

	cd([lastdir 'res/']);

	nmrpipe=['nmrpipe_dir/'];
	for rowi=1:nsample
	  sampdir=num2str(rowi);
	  oridir=rawfolders{rowi};
	  mkdir(['./' nmrpipe sampdir]);
	  mkdir([ './' nmrpipe sampdir '/script']);
      pause(1)
%       fprintf('/n %s /n',[datadir oridir '/*'])
%       fprintf('%s /n',[nmrpipe sampdir '/'])
%       fprintf('%s /n',[nmrpipe 'script'])
%       fprintf('%s /n',[nmrpipe sampdir '/script'])
	  copyfile([datadir oridir '/*'],[nmrpipe sampdir '/']);
      system(['rm -f -r ' lastdir 'res/' nmrpipe sampdir  '/SAND']);
	  copyfile([nmrpipe 'script'],[nmrpipe sampdir '/script']);
      system(['mv ' nmrpipe sampdir '/' 'test.ft' ' ' nmrpipe sampdir '/' 'raw.ft']);
	end
end

%%%-----change
% if (str2num(getenv('PREP')) == 0) % need to be changed as input data is now different
% 	cd(datadir);
% 	rawfiles=strsplit(strtrim(ls));
% 	nsample=size(rawfiles,2);
% 
% 	cd([lastdir '/res/']);
% 
% 	nmrpipe=['nmrpipe_dir/'];
% 	for rowi=1:nsample
% 	  sampdir=num2str(rowi);
% 	  orifile=rawfiles{rowi};
%       mkdir([nmrpipe sampdir]);
%       system(['rm -r ' nmrpipe sampdir]);
% 	  mkdir([nmrpipe sampdir]);
% 	  mkdir([nmrpipe sampdir '/script']);
% 	  copyfile([datadir orifile],[nmrpipe sampdir '/']);
% 	  system(['mv ' nmrpipe sampdir '/' orifile ' ' nmrpipe sampdir '/' 'raw.ft']);
% 	  copyfile([nmrpipe 'script'],[nmrpipe sampdir '/script']);
% 	end
% end
%%%----end change


% nmrpipe based preprocess
nmrpipe_process2('./nmrpipe_dir/',nsample,'prior',shelladd,str2num(getenv('PREP')));

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
  ppmrange=[str2num(getenv('LPPM')) str2num(getenv('UPPM'))];
  binrange=buckets.optimized.expandedBuckets;
  binrange=flip(binrange,2);
  binrange=binrange(binrange(:,1)<ppmrange(2)&binrange(:,2)>ppmrange(1),:);%limited ppm searching range
  % group every three bins
  binrange=binrange';
  binrange_vec=binrange(:);
  %combine each neighbor 3 range
  groupind=[];
  stepsize=2*str2num(getenv('GBIN'));
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

if (str2num(getenv('PREP')) == 1) 
	cd(datadir);
    nsample=size(rawfolders,2);   
    rawfolders=strsplit(strtrim(ls));
    for rowi=1:nsample
        cd(datadir)
        cd(rawfolders{rowi})
        system('rm -f -r SAND/last_preprocess');
        mkdir('SAND/last_preprocess');
 
    end	
    
	cd([lastdir 'res/']);
    nmrpipe=['nmrpipe_dir/'];
	for rowi=1:nsample
	  sampdir=num2str(rowi);
	  oridir=rawfolders{rowi};
	  copyfile([nmrpipe sampdir '/*'],[datadir oridir '/SAND/last_preprocess/']);
      %system(['mv ' nmrpipe sampdir '/' 'test.ft' ' ' nmrpipe sampdir '/' 'raw.ft']);
    end
    writetable([table((1:length(rawfolders))','VariableNames',"after"),cell2table(rawfolders', "VariableNames", "before")],[lastdir 'inputfolders.txt'],'Delimiter','\t');

end

