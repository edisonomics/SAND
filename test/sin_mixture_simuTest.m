function tests=sin_mixture_simuTest
  tests=functiontests(localfunctions);
end

function testsimplerun(testCase)
  load('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/code/testdata/sin_funct_exp.mat');
  para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
  timevec=1:0.1:10;
  resvec=sin_mixture_simu(para,timevec,nan,'real');
  boolflag=sum(resvec-resvec_store)<0.0001;
  verifyTrue(testCase,boolflag);
end

function testcomplex_real_equal(testCase)
  para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
  timevec=1:0.1:10;
  resvec=sin_mixture_simu(para,timevec,nan,'real');
  resvec_comp=sin_mixture_simu(para,timevec,nan,'complex');
  verifyEqual(testCase,imag(resvec_comp),resvec)
end

function testgausiannoise(testCase)
  para=[1 0 1 0; 2 0.3 2 0; 5 0.3 1 0];
  timevec=1:0.002:10;
  resvec_comp_noise=sin_mixture_simu(para,timevec,0.1,'complex');
  resvec_comp_nonoise=sin_mixture_simu(para,timevec,nan,'complex');
  sigvec=resvec_comp_noise-resvec_comp_nonoise;
  sigvec_real=real(sigvec);
  sigvec_imag=imag(sigvec);
  tolerance=0.01;
  realmean_0=abs(mean(sigvec_real))<tolerance;
  realstd=abs(std(sigvec_real)-0.1)<tolerance;
  imagmean_0=abs(mean(sigvec_imag))<tolerance;
  imagstd=abs(std(sigvec_imag)-0.1)<tolerance;
  verifyTrue(testCase,realmean_0&realstd&imagmean_0&imagstd);
end
