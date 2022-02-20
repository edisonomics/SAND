function [obj]=sum_sq_obj(para,realval,nsig,timevec,simutype)
% This function is the objective function for local optimization. It is used currently only in the local step in the optimization (because of the specific need of the format, and the global optimization objective function is coded inside of the function).
%
% Arguments:
%         para: vector. numeric. The parameter for simulating sin function. Parameters in order: f, lambda, A, phi. Parameters for each composing sin function are concatenated. Must be provided
%         realval: vector. numeric. the true value to compare with in the optimization. Must be provided.
%         nsig: integer. the number of signal in simulation. Used to reshape the para vector. Must be provided
%         timevec: vector. numeric. The time grid in the simulation. Must be provided
%         simutype: string. Either 'real' or 'complex'. Default 'real'. imag(X1)==X2 X1 is simulation with type=="complex" and X1 is with type=="real" (all other parameters the same).
% Return:
%         the objective function by comparing simulated spectra and real spectra
%
% Example:
%
% rng(1)
% para=[2 0.3 2 0];%f, lambda, A, phi
% timevec=1:0.0002:10;
% sdval=1;
% nsig=1;
% specsig=sin_mixture_simu(para,timevec,sdval,'complex');
% sse=sum_sq_obj(para,specsig,nsig,timevec,'complex');
% Test:
%   results = runtests('sum_sq_objTest.m')
% Yue Wu 07/14/2020

if ~exist('para','var')
  error('please provide input parameter vector');
end
if ~exist('realval','var')
  error('please provide compared vector');
end
if ~exist('nsig','var')
  error('please provide number of signals in the mixture');
end
if ~exist('timevec','var')
  error('please provide input time vector');
end
if ~exist('simutype','var')
  simutype='real';
end
% convert para from vec to array
nparadim=floor(length(para)/nsig);
pred_vec=sin_mixture_simu(reshape(para,[nparadim,nsig])',timevec,nan,'complex');
obj=sum(abs(realval-pred_vec).^2);
