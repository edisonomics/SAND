function tests=sum_sq_objTest
  tests=functiontests(localfunctions);
end

% this should original be tested with different sample size as sanity check. However, the standard error for sigma is surprisingly complex so ignored.
function testsd(testCase)
  rng(1)
  para=[2 0.3 2 0];%f, lambda, A, phi
  timevec=1:0.0002:10;
  sdval=1;
  nsig=1;
  specsig=sin_mixture_simu(para,timevec,sdval,'complex');
  sse=sum_sq_obj(para,specsig,nsig,timevec,'complex');
  equalflag=(sqrt(sse/length(timevec)/2)-sdval)/sdval<0.01;%/2 cause for complex data the sse is doubled
  verifyTrue(testCase,equalflag);
end
