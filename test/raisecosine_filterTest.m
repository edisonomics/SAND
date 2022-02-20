function tests=raisecosine_filterTest
  tests=functiontests(localfunctions);
end

% test rectangles
function test_rect(testCase)
  ppm=-1:0.01:11;
  ppmregion=[1.2 1.4];
  filt_beta=0;
  filvec=raisecosine_filter(ppm,ppmregion,filt_beta);
  rangein=matchPPMs(ppmregion,ppm);
  in_ind=rangein(1):rangein(2);
  ppmmask=false([1,size(ppm,2)]);
  ppmmask(in_ind)=true;
  plateau_flag=all(filvec(ppmmask)==1);
  baseline_flag=all(filvec(~ppmmask)==0);
  verifyTrue(testCase,plateau_flag&baseline_flag);
end

% test softbounds
% 1 at (1-beta)/2T 0.5 at 1/2T 0 at (1+beta)/2T
function test_bounds(testCase)
  ppm=-1:0.01:11;
  ppmregion=[2 4];
  filt_beta=0.2;
  filvec=raisecosine_filter(ppm,ppmregion,filt_beta);
  midppm=mean(ppmregion);
  halfb=(ppmregion(2)-ppmregion(1))/2;
  % 1 at plateau range
  plateau_b_ind=matchPPMs(midppm-halfb*(1-filt_beta),ppm);
  flag_plateau=(filvec(plateau_b_ind)==1);
  % 0.5 at 1/2T
  half_ind=matchPPMs(midppm-halfb,ppm);
  flag_half=abs(filvec(half_ind)-0.5)<0.000001;
  % 0 at (1+beta)/2T
  tail_ind=matchPPMs(midppm-halfb*(1+filt_beta),ppm);
  flag_tail=abs(filvec(tail_ind)-0)<0.000001;
  verifyTrue(testCase,flag_plateau&flag_half&flag_tail);
end
