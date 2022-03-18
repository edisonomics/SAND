#!/bin/csh

bruk2pipe -verb -in ./fid \
  -bad 0.0 -ext -aswap -AMX -decim 1664 -dspfvs 21 -grpdly 76  \
  -xN             65536  \
  -xT             32691  \
  -xMODE            DQD  \
  -xSW        12019.231  \
  -xOBS         600.133  \
  -xCAR           4.754  \
  -xLAB              1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 1.22070e-01 \
  -out ./test2.fid -ov

sleep 5
