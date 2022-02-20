function [filt_vec]=raisecosine_filter(ppm,ppmregion,filt_beta)
% this function will produce a filter vector for selected ppm region, following raise cos filter
% The filter will have plateau in the middle and have soft boundaries
% More information here https://en.wikipedia.org/wiki/Raised-cosine_filter
% Arguments:
%         ppm: vector. numeric. the ppm vector in NMR spectra. Must be provided
%         ppmregion: vector. numeric. two element. the region to have the spectral filter (in terms of filt_beta=0). Must be provided
%         filt_beta: numeric. beta parameter in the filter equation. Default 0.
% Return:
%         The filter vector that can be multiplied with the original spectral signal
% Example:
%
% ppm=-1:0.01:11;
% ppmregion=[1.2 1.4];
% filt_beta=0;
% filvec=raisecosine_filter(ppm,ppmregion,filt_beta);
% plot(ppm,filvec);
% ylim([-1,2]);
% ppmregion=[2 4];
% filt_beta=0.2;
% filvec=raisecosine_filter(ppm,ppmregion,filt_beta);
% plot(ppm,filvec);
% ylim([-1,2]);
%
% Test:
%   results = runtests('raisecosine_filterTest.m')
% Yue Wu 09/14/2020

if ~exist('ppm','var')
  error('please provide input ppm vector');
end
if ~exist('ppmregion','var')
  error('please provide input selected filter ppm range');
end
if ~exist('filt_beta','var')
  filt_beta=0;%hard boundaries as in rectangles
end
midppm=mean(ppmregion);
filtf_step=(ppmregion(2)-ppmregion(1))/2;%1/2T
lenvec=size(ppm,2);
filt_vec=zeros([1,lenvec]);
% plateau
del_plateau=filtf_step*(1-filt_beta);
rang_plateau=[midppm-del_plateau,midppm+del_plateau];
oneinds=matchPPMs(rang_plateau,ppm);
filt_vec(oneinds(1):oneinds(2))=1;
% shoulders
filtT=1/(filtf_step*2);
f_shoulder_rela_rang=[del_plateau,(filtf_step*(1+filt_beta))];
shoulderinds=matchPPMs(f_shoulder_rela_rang+midppm,ppm);
if shoulderinds(1)~=shoulderinds(2)
  f_shoulder_rela_vec=ppm(shoulderinds(1):shoulderinds(2))-midppm-del_plateau;
  coscurv=(1+cos(pi*filtT/filt_beta*f_shoulder_rela_vec))/2;
  %% right side
  lowrange=[midppm-f_shoulder_rela_rang];
  rngind=matchPPMs(lowrange,ppm);
  filt_vec(rngind(2):rngind(1))=flip(coscurv);
  %% left side
  highrange=[midppm+f_shoulder_rela_rang];
  rngind=matchPPMs(highrange,ppm);
  filt_vec(rngind(1):rngind(2))=coscurv;
end
