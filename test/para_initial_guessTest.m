function tests=para_initguessTest
  tests=functiontests(localfunctions);
end

function testexample(testCase)
  para=[1 0.3 1 0; 2 0.3 2 0; 5 0.3 1 0];
  timevec=1:0.01:10;
  residual_sig=sin_mixture_simu(para,timevec,0.1,'complex');
  fseq=0.5:0.5:6;
  lambda_default=0.3;
  rela_range=0.1;
  [para_initguess para_iniguess_rang]=para_initial_guess(residual_sig,timevec,fseq,lambda_default,rela_range);
  flag_f=para_initguess(1)==2.0;
  meanval=mean(para_iniguess_rang,2);
  down_ratio=abs(para_iniguess_rang(:,1)-meanval)./meanval;
  up_ratio=abs(para_iniguess_rang(:,2)-meanval)./meanval;
  flag_relat=all(abs(down_ratio(1:3)-rela_range)<0.00001 & abs(up_ratio(1:3)-rela_range)<0.00001);
  flag_unchange=abs(para_iniguess_rang(4,1)-para_iniguess_rang(4,2))<0.00001;
  verifyTrue(testCase,flag_f&flag_relat&flag_unchange);
end
