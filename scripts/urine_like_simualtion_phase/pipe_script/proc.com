#!/bin/csh

cd ..
# read the saved fid from Matlab saved files
txt2pipe.tcl -x -complex -time -in temp.fid.txt -inHdr org.fid > temp.fid

nmrPipe -in temp.fid \
| nmrPipe  -fn ZF -auto \
| nmrPipe  -fn FT -auto -verb \
-out test.ft1 -ov -xi

# change ft header
sethdr test.ft1 -title spec
cd ./script/
