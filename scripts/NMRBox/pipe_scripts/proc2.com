#!/bin/csh
# convert the fid into nmrpipe format
# NMRBox default NMRPipe function for already preprocessed data (need to be changed by USER based on data): OS 07-04-2022
# fd 1/13/2023 replace (ref1D, interpNMR, sethdr -title) with copyHdr.tcl

cd ..

echo "Executing SAND NMRPipe script proc2.com (`pwd`) ..."

#bruker -AUTO > conv.out
#chmod a+rx ./fid.com
#fid.com > conv.out

nmrPipe -in raw.ft \
| nmrPipe -fn HT \
| nmrPipe -fn FT -inv \
| nmrPipe -fn ZF -inv \
 -out test.fid -ov

#
# fd and omid 12/21/2022: 
#  Add -xC 1.0 to correct for first point scaling.
#  Save scale factor (was -scaleTo 1000.0)

# 1/31/2024 add -xELB 0.0 fd lp and zh
# basicFT1.com -in test.fid -out test.ft1 -xC1 1.0
basicFT1.com -in test.fid -out test.ft1 -xC1 1.0 -xELB 0.0

set maxVal      = (`specStat.com -in test.ft1 -stat vMax -brief`)
set scaleFactor = (`MATH "div( 1000.0, $maxVal )"`)

nmrPipe -in test.ft1 -out test.ft1 -fn MULT -c $scaleFactor -inPlace

echo $scaleFactor > scale.txt

#
# fd 1/13/2023
#ref1D.tcl -in test.ft1 -tab ref.tab -x1 0.05ppm -xn -0.05ppm -xRef 0.0 -max
#interpNMR -ref ../1/test.ft1 -in test.ft1 -out test.ft1

copyHdr.tcl -in ../1/test.ft1 -out test.ft1

# 1/31/2024 add -xC1 1.0 -xELB 0.0 fd lp and zh
# basicFT1.com -in test.fid -out test2.ft1
basicFT1.com -in test.fid -out test2.ft1 -xC1 1.0 -xELB 0.0



set maxVal      = (`specStat.com -in test2.ft1 -stat vMax -brief`)
set scaleFactor = (`MATH "div( 1000.0, $maxVal )"`)

nmrPipe -in test2.ft1 -out test2.ft1 -fn MULT -c $scaleFactor -inPlace

echo $scaleFactor > scale2.txt

interpNMR -ref ../1/test.ft1 -in test2.ft1 -out test2.ft1
# rescale.com -in test.ft1 -out test.ft1 -x1 2.8ppm -xn 1.8ppm

# fd 1/13/2023
# change ft header
# sethdr test.ft1 -title spec

cd ./script/
