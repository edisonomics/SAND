#!/bin/csh

# clean up the folder
genericClean.com -fList nmr.com nmr2.com
cd ../
genericClean.com -dList mask mask_fid mask_ft -eSkip .fid
rm -r  *_matlab
cd ./script
