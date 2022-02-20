#!/bin/csh

# clean up the folder
cd ../
genericClean.com -dList mask mask_fid mask_ft -eSkip .fid
rm -r  *_matlab
cd ./script
