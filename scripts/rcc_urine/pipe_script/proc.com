#!/bin/csh

cd ..
# read the saved ft txt from Matlab saved files
txt2pipe.tcl -x -real -freq -in temp.ft.txt -xN 68822 -xT 34411 -xMODE Real -xSW 6301.3 -xOBS 600.133  -xCAR 4.754 -xFT Freq -xLAB H1 > test.ft1

cd ./script/
