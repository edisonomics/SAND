#!/bin/csh
# convert the fid into nmrpipe format
# NMRBox default NMRPipe function (need to be changed by USER based on data): OS 07-18-2022
# fd 1/13/2023 replace (ref1D, interpNMR, sethdr -title) with copyHdr.tcl

cd ..

echo "Executing SAND NMRPipe script proc.com (`pwd`) ..."

bruker -AUTO > conv.out
chmod a+rx ./fid.com
fid.com > conv.out
#
basicFT1.com \
 -in test.fid -scaleTo 1000.0 -xELB 0.3 \
 -xBASEARG POLY,auto,ord=0,window=2% -xP0 Auto -xP1 0 \
 -apOrd 0 -apArgs apx1=-1ppm,apxn=-3ppm,apMode=rms,apxP0Step=0.1  -out test.ft1

# a process without line broadening
basicFT1.com \
 -in test.fid -scaleTo 1000.0 -xELB 0 \
 -xBASEARG POLY,auto,ord=0,window=2% -xP0 Auto -xP1 0 \
 -apOrd 0 -apArgs apx1=-1.0ppm,apxn=-3.0ppm,apMode=rms,apxP0Step=0.1 -out test2.ft1

#
# fd 1/13/2023 (rescale step was already removed).
# ref1D.tcl -in test2.ft1 -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
# interpNMR -ref ../1/test.ft1 -in test2.ft1 -out test2.ft1
# rescale.com -in test.ft1 -out test.ft1 -x1 2.8ppm -xn 1.8ppm
# ref1D.tcl -in test.ft1 -tab ref.tab -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
# interpNMR -ref ../1/test.ft1 -in test.ft1 -out test.ft1

copyHdr.tcl -in ../1/test.ft1 -out test.ft1
copyHdr.tcl -in ../1/test.ft1 -out test2.ft1

# fd 1/13/2023
# change ft header
# sethdr test.ft1 -title spec

cd ./script/
