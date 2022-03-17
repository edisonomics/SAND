#!/bin/csh

var2pipe -verb -in ../fid \
 -noaswap  \
  -xN             32768  \
  -xT             16384  \
  -xMODE        Complex  \
  -xSW         9000.900  \
  -xOBS         499.409  \
  -xCAR           4.773  \
  -xLAB              H1  \
  -ndim               1  \
| nmrPipe -fn MULT -c 1.25000e+02 \
  -out ../test.fid -ov
