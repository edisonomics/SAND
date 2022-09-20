#!/bin/csh

# clean up the folder
# NMRBox default NMRPipe function (don't need to be changed): OS 6-20-2022
genericClean.com -fList nmr.com nmr2.com
cd ../
genericClean.com -dList mask mask_fid mask_ft -eSkip .fid
rm -f -r  *_matlab
cd ./script
