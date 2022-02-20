function tests=para_initial_guess_fftTest
  tests=functiontests(localfunctions);
end

% the test depends on local files
function testsimplerun(testCase)
  cd('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/test');
  datadir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/band_generate_run_locpara_precis/';
  testfid='ori_matlab/test_trans_ift.txt';
  testfidhdr='test_trans.fid';
  testft='test.ft1';
  % load ft
  [value axes]=read_nmrp([datadir testft]);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1;
  fttab=table(ppm,value);
  [maxint,maxind]=max(fttab{:,2});
  % load fid
  fidtab=readtable([datadir testfid],'Format','%f%f%f');
  fidtab.Properties.VariableNames={'time','real','imag'};
  residual_sig=fidtab{:,2}+i*fidtab{:,3};
  hdrpath=[datadir testfidhdr];
  conv_f=[axes(8,1) axes(14,1)];
  rela_range=0.1;
  [para_initguess para_iniguess_rang]=para_initial_guess_fft(residual_sig.',hdrpath,conv_f,rela_range);
  expf=(fttab{maxind,1}-conv_f(1))*conv_f(2);
  f_same_flag=abs(expf-para_initguess(1))<0.1;
  % A_same_flag=(maxint-para_initguess(3))<1;
  verifyTrue(testCase,f_same_flag);
end

% the test depends on local files
function testsimplerun2(testCase)
  cd('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/test');
  datadir='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spec_deconv_time_domain/result/band_generate_run_locpara_precis/';
  testfid='mask_fid_matlab/test251.fid.txt';
  testfidhdr='mask_fid/test251.fid';
  testft='mask_ft/test251.ft1';
  % load ft
  [value axes]=read_nmrp([datadir testft]);
  ppm=inc2ppm(axes);
  ppm=ppm.ppm1;
  fttab=table(ppm,value);
  [maxint,maxind]=max(fttab{:,2});
  % load fid
  fidtab=readtable([datadir testfid],'Format','%f%f%f');
  fidtab.Properties.VariableNames={'time','real','imag'};
  residual_sig=fidtab{:,2}+i*fidtab{:,3};
  hdrpath=[datadir testfidhdr];
  conv_f=[axes(8,1) axes(14,1)];
  rela_range=0.1;
  [para_initguess para_iniguess_rang]=para_initial_guess_fft(residual_sig.',hdrpath,conv_f,rela_range);
  expf=(fttab{maxind,1}-conv_f(1))*conv_f(2);
  f_same_flag=abs(expf-para_initguess(1))<0.1;
  % A_same_flag=(maxint-para_initguess(3))<1;
  verifyTrue(testCase,f_same_flag);
end
