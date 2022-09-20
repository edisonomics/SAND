function [fttab]=ft_pipe3(fidtab,hd_file,dirname,addstr)
% fourier transformation as implemeted in nmrpipe for specific usage if NMR data
%
% Arguments:
%         fidtab: table. numeric. two colume: time, fid (complex number). Must be provided
%         hd_file: str. the location of the original header for the fid file. Must be provided
%         dirname: name of the temp folder. Must be provided
%         addstr: str. string to addon before nmrpipe running in the shell. mostly used for wired shell environments
% Return:
%         transformed ft signal table (ppm,intensity). The output file and shell script will be created locally
% Examples:
%
% cd('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/test');
% datadir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/band_generate_run_locpara_precis/';
% testfid='ori_matlab/test_trans_ift.txt';
% testfidhdr='test_trans.fid';
% testft='test.ft1';
% [value axes]=read_nmrp([datadir testft]);
% ppm=inc2ppm(axes);
% ppm=ppm.ppm1;
% fttab=table(ppm,value);
% fidtab=readtable([datadir testfid],'Format','%f%f%f');
% fidtab.Properties.VariableNames={'time','real','imag'};
% fttabnew=ft_pipe(fidtab,[datadir testfidhdr],'temp');
% plot(fttab{:,2}-fttabnew{:,2});
% plotr(fttab{:,1},fttab{:,2});
% hold on;
% plotr(fttabnew{:,1},fttabnew{:,2});
%
% Tests:
%
% results = runtests('ft_pipeTest.m')
% Yue Wu 11/02/2020
% OS 18/06/2022 (revised)

if ~exist('fidtab','var')
  error('please provide input fid array');
end
if ~exist('hd_file','var')
  error('please provide the input header');
end
if ~exist('dirname','var')
  error('please provide the folder name');
end
if ~exist('addstr','var')
  addstr='';
end
% output fid
% arrayoutput(:,1)=fidtab(:,1);
mkdir(dirname);
cd(dirname);
ind=1:size(fidtab,1);
realpart=fidtab{:,2};
imagpart=fidtab{:,3};
taboutput=table(ind',realpart,imagpart);
writetable(taboutput,'temp.fid.txt','Delimiter',' ','WriteVariableNames',false);
copyfile(hd_file,['./org.fid']);
% check current shell of the system
shell_line="#!/bin/csh";
[stat out]=system('echo $0');
if contains(out,'tcsh')
  shell_line="#!/bin/tcsh";
end
% create ft_pipe.sh file
shellcommand=[shell_line,
              addstr,
              ['txt2pipe.tcl -x -complex -time -in temp.fid.txt -inHdr org.fid > temp.fid'],
              "nmrPipe -in temp.fid \\",
              "| nmrPipe  -fn ZF -auto \\",
              "| nmrPipe  -fn FT -auto -verb \\",
              "  -out temp.ft -ov -xi"];
shellcommand_comb=join(shellcommand,'\n');
FID=fopen('ft_pipe.sh','w');
fprintf(FID,shellcommand_comb);
fprintf(FID,'\n');
fclose(FID);


%fprintf(fopen('ft_pipe.sh','wt'),fopen('ft_pipe.sh'));
% run ft_pipe.sh
system('chmod +x ./ft_pipe.sh');
system('./ft_pipe.sh');
% load data
[value axes]=read_nmrp('temp.ft');
ppm=inc2ppm(axes);
ppm=ppm.ppm1;
fttab=table(ppm,value);
cd('../');
% if size(fttab,1)~=0
%   rmdir(dirname,'s');
% end
