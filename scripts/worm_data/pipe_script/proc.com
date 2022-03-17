#!/bin/csh

cd ..
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
cd ./script/
