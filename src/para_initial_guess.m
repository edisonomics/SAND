function [para_initguess para_iniguess_rang]=para_initial_guess(residual_sig,timevec,fseq,lambda_default,rela_range)
% Produce initial guess for paramters based on the residual signal (original_signal - current_signal). Parameter initial guesses for one new signal will be returend
%
% Arguments:
%         residual_sig: vector. numeric. The time domain signal residual (subsampled). dimension [1 #time_point]. Must be provided.
%         timevec: vector. numeric. The time grid for the time domain signal. dimension [1 #time_point]. Must be provided.
%         fseq: vector. numeric. The frequncy sequence to sweep through. choice of value depends on the range of frequency. Must be provided.
%         lambda_default: struct. The default lambda value. Must be provided.
%         rela_range: numeric. The relative uncerntainty for paramter searching range. The calculated range will be [(1-rela_range)*initial_guess (1+rela_range)*initial_guess]. Must be provided
% Return:
%         para_initguess: vector. numeric. The initial guess parameter the next add-on model. each element is a different parameter:  f, lambda, A, phi. Initial guess of phi is always 0 and so the range will be 0 as well.
%         para_iniguess_rang: array. numeric. The parameter range for the next add-on model. each row is a different parameter:  f, lambda, A, phi. npara*2 npara=4
% Examples:
%
% para=[1 0.3 1 0; 2 0.3 2 0; 5 0.3 1 0];
% timevec=1:0.01:10;
% residual_sig=sin_mixture_simu(para,timevec,0.1,'complex');
% fseq=0.5:0.5:6;
% lambda_default=0.3;
% rela_range=0.1;
% [para_initguess para_iniguess_rang]=para_initial_guess(residual_sig,timevec,fseq,lambda_default,rela_range)
%
% Test:
% results = runtests('para_initial_guessTest.m')
% Yue Wu 06/23/2020
if ~exist('residual_sig','var')
  error('please provide input time domain signal residual');
end
if ~exist('timevec','var')
  error('please provide input time grid');
end
if ~exist('fseq','var')
  error('please provide input frequncy grid');
end
if ~exist('lambda_default','var')
  error('please provide input default lambda value');
end
if ~exist('rela_range','var')
  error('please provide relative range for parameters');
end
Avec=[];
for fele=fseq
  para_simu=[fele lambda_default 1 0];%f, lambda, A, phi. The simulation is based on default value for lambda and phi, so will be affected if these two cannot be correctly guessed.
  s_vec=sin_mixture_simu(para_simu,timevec,nan,'complex');
  residual_flat=[real(residual_sig) imag(residual_sig)];
  s_flat=[real(s_vec) imag(s_vec)];
  Aele=abs(sum(residual_flat.*s_flat)/sum(s_flat.*s_flat));
  Avec=[Avec Aele];
end
ini_ind=find(Avec==max(Avec));
para_initguess=[fseq(ini_ind) lambda_default max(Avec) 0];
para_iniguess_rang=[];
for i=1:length(para_initguess)
  para_iniguess_rang=[para_iniguess_rang; para_initguess(i)*[1-rela_range 1+rela_range]];
end
