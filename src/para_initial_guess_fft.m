function [para_initguess para_iniguess_rang]=para_initial_guess_fft(residual_sig,hdrpath,conv_f,rela_range)
% Produce initial guess for paramters based on the residual signal (original_signal - current_signal). Parameter initial guesses for one new signal will be returend. This function will use fft in nmrpipe to find the maximal A
%
% Arguments:
%         residual_sig: vector. numeric. The time domain signal residual (subsampled). dimension [1 #time_point]. Must be provided.
%         hdrpath: string. The location of the original hdr file for fid. Must be provided.
%         conv_f: vector. numeric. The conversion factor from ppm to frequncy. Must be provided.
%         rela_range: numeric. The relative uncerntainty for paramter searching range. The calculated range will be [(1-rela_range)*initial_guess (1+rela_range)*initial_guess]. Must be provided
% Return:
%         para_initguess: vector. numeric. The initial guess parameter the next add-on model. each element is a different parameter:  f, lambda, A, phi. Initial guess of phi is always 0 and so the range will be 0 as well.
%         para_iniguess_rang: array. numeric. The parameter range for the next add-on model. each row is a different parameter:  f, lambda, A, phi. npara*2 npara=4
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
% [maxint,maxind]=max(fttab{:,2});
% fidtab=readtable([datadir testfid],'Format','%f%f%f');
% fidtab.Properties.VariableNames={'time','real','imag'};
% residual_sig=fidtab{:,2}+i*fidtab{:,3};
% hdrpath=[datadir testfidhdr];
% conv_f=[axes(8,1) axes(14,1)];
% rela_range=0.1;
% [para_initguess para_iniguess_rang]=para_initial_guess_fft(residual_sig.',hdrpath,conv_f,rela_range);
% expf=(fttab{maxind,1}-conv_f(1))*conv_f(2);
%
% Test:
% results = runtests('para_initial_guess_fftTest.m')
% Yue Wu 06/23/2020
if ~exist('residual_sig','var')
  error('please provide input time domain signal residual');
end
if ~exist('hdrpath','var')
  error('please provide input header file path');
end
if ~exist('conv_f','var')
  error('please provide input conversion factor');
end
if ~exist('rela_range','var')
  error('please provide relative range for parameters');
end
%
residual_sig=residual_sig.';
spec_new=ft_pipe(table([1:length(residual_sig)]',real(residual_sig),imag(residual_sig)),hdrpath,'temp');
inthere=spec_new{:,2};
ppmhere=spec_new{:,1};
[Aguss,Aind]=max(inthere);
ppmguess=ppmhere(Aind,1);
fguess=(ppmguess-conv_f(1))*conv_f(2);
halfh_ind=find((inthere-Aguss/2)<0.00001);
[~,sortind]=sort(abs(halfh_ind-Aind));
lambda_ppm=abs(ppmhere(halfh_ind(sortind(1)))-ppmhere(Aind))*2;
lambda_default=lambda_ppm*conv_f(2);
para_initguess=[fguess lambda_default Aguss 0];
para_iniguess_rang=[];
for i=1:length(para_initguess)
  para_iniguess_rang=[para_iniguess_rang; para_initguess(i)*[1-rela_range 1+rela_range]];
end
