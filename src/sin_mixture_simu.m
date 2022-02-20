function [resvec]=sin_mixture_simu(para,timevec,sigma,type)
% simulating decay sin function mixture based on provided parameters
%
% Arguments:
%         para: array. numeric. The parameter for simulating sin function. each row is a different component, and each complumn is for one parameter. Parameters in order: f, lambda, A, phi. The number of rows will be the number of composing sin functions. Must be provided
%         timevec: vector. numeric. The time grid in the simulation. Must be provided
%         sigma: numeric. The normal sampled noise for the mixture signal. Must be provided
%         type: string. Either 'real' or 'complex'. Default 'real'. imag(X1)==X2 X1 is simulation with type=="complex" and X1 is with type=="real" (all other parameters the same).
% Return:
%         resvec the simualted sin function mixture
% Examples:
%
% para=[1 0 1 0];
% timevec=1:0.01:5;
% resvec=sin_mixture_simu(para,timevec,nan,'real');
% plot(timevec,resvec);
%
% para=[1 0 1 0; 2 0.1 2 0];
% timevec=1:0.01:5;
% resvec=sin_mixture_simu(para,timevec,0.1,'real');
% plot(timevec,resvec);
%
% para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
% timevec=1:0.01:10;
% resvec=sin_mixture_simu(para,timevec,0.1,'real');
% plot(timevec,resvec);
%
% rng(1)
% para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
% timevec=1:0.01:10;
% resvec=sin_mixture_simu(para,timevec,nan,'real');
% plot(timevec,resvec);
%
% para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
% timevec=1:0.01:10;
% resvec=sin_mixture_simu(para,timevec,nan,'complex');
% figure();
% plot(timevec,real(resvec));
% figure();
% plot(timevec,imag(resvec));
%
% para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
% timevec=1:0.01:10;
% resvec=sin_mixture_simu(para,timevec,0.1,'complex');
% figure();
% plot(timevec,real(resvec));
% figure();
% plot(timevec,imag(resvec));
%
% Tests:
% results = runtests('sin_mixture_simuTest.m')
% Yue Wu 06/19/2020

if ~exist('para','var')
  error('please provide input parameter array');
end
if ~exist('timevec','var')
  error('please provide the time grid');
end
if ~exist('sigma','var')
  error('please provide input sigma for the noise distribution');
end
if ~exist('type','var')
  type='real';
end

resvec=zeros(size(timevec));
for rowi=1:size(para,1)
  paravec=para(rowi,:);
  if strcmp(type,'real')
    resvec=resvec+(paravec(3).*sin(2*pi*paravec(1)*timevec+paravec(4))).*exp(-paravec(2)*timevec); %f, lambda, A, phi
  elseif strcmp(type,'complex')
    resvec=resvec+paravec(3).*exp(1i*(2*pi*paravec(1)*timevec+paravec(4))-paravec(2)*timevec);
  end
end
if ~isnan(sigma)
  resvec=resvec+normrnd(0,sigma,size(resvec));
  if strcmp(type,'complex')
    resvec=resvec+1i*normrnd(0,sigma,size(resvec));%%gaussian noise in both real and imaginary part
  end
end
