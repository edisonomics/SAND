function tests=ft_pipeTest
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
  % load fid
  fidtab=readtable([datadir testfid],'Format','%f%f%f');
  fidtab.Properties.VariableNames={'time','real','imag'};
  fttabnew=ft_pipe(fidtab,[datadir testfidhdr],'temp');
  % plot(fttab{:,2}-fttabnew{:,2});
  % plotr(fttab{:,1},fttab{:,2});
  % hold on;
  % plotr(fttabnew{:,1},fttabnew{:,2});
  boolflag=abs(sum(fttab{:,2}-fttabnew{:,2}))<0.1;
  verifyTrue(testCase,boolflag);
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
  % load fid
  fidtab=readtable([datadir testfid],'Format','%f%f%f');
  fidtab.Properties.VariableNames={'time','real','imag'};
  fttabnew=ft_pipe(fidtab,[datadir testfidhdr],'temp');
  % plot(fttab{:,2}-fttabnew{:,2});
  % plotr(fttab{:,1},fttab{:,2});
  % hold on;
  % plotr(fttabnew{:,1},fttabnew{:,2});
  boolflag=abs(sum(fttab{:,2}-fttabnew{:,2}))<0.3;
  verifyTrue(testCase,boolflag);
end
