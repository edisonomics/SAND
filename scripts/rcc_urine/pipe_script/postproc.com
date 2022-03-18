#!/bin/csh

# transform nmrpipe files for matlab input
# transform the fid, ft files to text table files

cd ../
rm -r mask_matlab mask_fid_matlab mask_ft_matlab sum_matlab ori_matlab
mkdir mask_matlab mask_fid_matlab mask_ft_matlab sum_matlab ori_matlab

cd ./mask_fid
foreach thefile (`ls *`)
  echo $thefile
  pipe2txt.tcl -index sec $thefile > ../mask_fid_matlab/${thefile}.txt &
end

cd ../
# summed spectra
pipe2txt.tcl -index PPM mask_sum.ft1 > ./sum_matlab/mask_sum.txt
pipe2txt.tcl -index sec mask_sum_ift.fid > ./sum_matlab/mask_sum_ift.txt
pipe2txt.tcl -index PPM mask_sum_ft.ft1 > ./sum_matlab/mask_sum_ft.txt

# original spectra
nmrPipe -in test.ft1 \
| nmrPipe -fn HT \
| nmrPipe -fn FT -inv \
| nmrPipe -fn ZF -inv \
 -out test_trans.fid -ov

pipe2txt.tcl -index PPM -fmt %e test.ft1 > ./ori_matlab/test_ft.txt
pipe2txt.tcl -index sec -fmt %e test_trans.fid > ./ori_matlab/test_trans_ift.txt

./script/fid.com
pipe2txt.tcl -index sec -fmt %e test2.fid > ./temp.fid.txt #without line broadening

cd ./script
