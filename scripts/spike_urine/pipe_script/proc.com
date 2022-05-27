#!/bin/csh

cd ..
# convert the fid into nmrpipe format
bruker -AUTO > conv.out
chmod a+rx ./fid.com
fid.com > conv.out
#
basicFT1.com -in test.fid -out test.ft1 -scaleTo 1000.0 -xELB 0.3 \
             -xP0 Auto -xP1 0 -xBASEARG POLY,auto,ord=0,window=2% \
             -apOrd 0 -apArgs apx1=12.0ppm,apxn=8.0ppm,apMode=rms,apxP0Step=0.1

ref1D.tcl -in test.ft1 -tab ref.tab -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
interpNMR -ref ../1/test.ft1 -in test.ft1 -out test.ft1
# a process without line broadening
basicFT1.com -in test.fid -out test2.ft1 -scaleTo 1000.0 -xELB 0 \
             -xP0 Auto -xP1 0 -xBASEARG POLY,auto,ord=0,window=2% \
             -apOrd 0 -apArgs apx1=12.0ppm,apxn=8.0ppm,apMode=rms,apxP0Step=0.1

ref1D.tcl -in test2.ft1 -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
interpNMR -ref ../1/test.ft1 -in test2.ft1 -out test2.ft1
# rescale.com -in test.ft1 -out test.ft1 -x1 2.8ppm -xn 1.8ppm

# change ft header
sethdr test.ft1 -title spec
cd ./script/
