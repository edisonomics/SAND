#!/bin/csh
# convert the fid into nmrpipe format
# NMRBox default NMRPipe function for already preprocessed data (need to be changed by USER based on data): OS 07-04-2022

cd ..

#bruker -AUTO > conv.out
#chmod a+rx ./fid.com
#fid.com > conv.out
#
nmrPipe -in raw.ft \
| nmrPipe -fn HT \
| nmrPipe -fn FT -inv \
| nmrPipe -fn ZF -inv \
 -out test.fid -ov

basicFT1.com -in test.fid -out test.ft1 -scaleTo 1000.0
ref1D.tcl -in test.ft1 -tab ref.tab -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
interpNMR -ref ../1/test.ft1 -in test.ft1 -out test.ft1

basicFT1.com -in test.fid -out test2.ft1 -scaleTo 1000.0 
ref1D.tcl -in test2.ft1 -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
interpNMR -ref ../1/test.ft1 -in test2.ft1 -out test2.ft1
# rescale.com -in test.ft1 -out test.ft1 -x1 2.8ppm -xn 1.8ppm

# change ft header
sethdr test.ft1 -title spec
cd ./script/
