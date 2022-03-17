#!/bin/csh

genericClean.com -fList nmr.com nmr2.com
cd ..
genericClean.com -dList mask mask_fid mask_ft
rm -r  *_matlab
cd ./script
