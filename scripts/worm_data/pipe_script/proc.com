#!/bin/csh
# generate ft script (in a simple way)
cd ../
#
# nmrPipe -in test.fid -fn SIGN -alt | nmrPipe -fn SIGN -i -out test.fid -inPlace
# nmrPipe -in test.fid -out test.fid -inPlace -fn MULT -xn 1 -c 2.0
#
basicFT1.com -in test.fid -out test.ft1 -scaleTo 1000.0 \
             -xELB 1.5 -list > ./script/nmr.com
chmod a+rx ./script/nmr.com
./script/nmr.com
# a process without line broadening
basicFT1.com -in test.fid -out test2.ft1 -scaleTo 1000.0 \
             -xELB 0 -list > ./script/nmr2.com
chmod a+rx ./script/nmr2.com
./script/nmr2.com
# change ft header
sethdr test.ft1 -title spec
cd ./script
