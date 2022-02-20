function [returnstate]=nmrpipe_process(pathfold,nfolder,stepflag)
% Run the different process scripts of NMRPipe
%
% Arguments:
%         pathfold: string. path to the folder to run process. Must be provided
%         nfolder: numeric. the number of folder to process. Must be provided
%         stepflag: string. the process stage. 'prior': the first step preprocessing. 'band': band filtering and other process. Must be provided
% Return:
%         0 if succeed. The files in the folder will be processed
% Examples:
%
% Tests:
%
% Yue Wu 02/20/2021

if ~exist('pathfold','var')
  error('please provide the directory path');
end
if ~exist('nfolder','var')
  error('please provide the number of folders to process');
end
if ~exist('stepflag','var')
  error('please provide process stage');
end
%
currdir=pwd;
cd(pathfold);
% check current shell of the system
shell_line="#!/bin/csh";
[stat out]=system('echo $0');
if contains(out,'tcsh')
  shell_line="#!/bin/tcsh";
end
if strcmp(stepflag,'prior')
  % create preprocess.sh file
  shellrun='preprocess.sh';
  shellcommand=[shell_line,
                'foreach thedir ( `seq 1 ' num2str(nfolder) '` )',
                '   echo $thedir',
                '   cd $thedir/script',
                '   ./preproc.com',
                '   cd ../../',
                'end'];
elseif strcmp(stepflag,'band')
  % create bandprocess.sh file
  shellrun='bandprocess.sh';
  shellcommand=[shell_line,
                'foreach thedir ( `seq 1 ' num2str(nfolder) '` )',
                '   echo $thedir',
                '   cd $thedir/script',
                '   ./band_postprepro.com',
                '   cd ../../',
                'end'];
end
shellcommand_comb=join(shellcommand,'\n');
FID=fopen(shellrun,'wt');
fprintf(FID,shellcommand_comb);
fclose(FID);
% run the shell
system(['chmod +x ./',shellrun]);
system(['./',shellrun]);
%
cd(currdir);
returnstate=0;
