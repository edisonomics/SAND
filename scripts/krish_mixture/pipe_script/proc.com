#!/bin/csh

# Conversion and processing of 1D data.
# This agilent-format raw data already has phase corrections recorded, so we can extract them.
cd ..
set vjArgs = (vj2pipe.com `vj2pipe.com -noverb -info |& tr '\\' ' '`)
set p0     = (`getArgD $vjArgs -xP0 0.0`)
set p1     = (`getArgD $vjArgs -xP1 0.0`)

# Create and execute a conversion script:
cd script
chmod a+rx fid.com
./fid.com
cd ..

# Process the data using the extracted phase corrections.
# Recalibrate ppm using vnmrj parameters rfl and rfp.
# The ppm calculation here assumes full spectral width (i.e. no extracted region).
basicFT1.com -in test.fid -out test.ft1 -xELB 1.5 -xP0 $p0 -xP1 $p1 -xBASEARG POLY,auto,ord=0 -scaleTo 1000.0 -list > ./script/nmr.com
chmod a+rx ./script/nmr.com
./script/nmr.com
sethdr test.ft1 -title MEASURED
# with no line broadening
basicFT1.com -in test.fid -out test2.ft1 -xELB 0 -xP0 $p0 -xP1 $p1 -xBASEARG POLY,auto,ord=0 -scaleTo 1000.0 -list > ./script/nmr2.com
chmod a+rx ./script/nmr2.com
./script/nmr2.com
# recalculate the ppm
set sw     = (`getParm -in test.ft1 -parm NDSW   -dim CUR_XDIM`)
set size   = (`getParm -in test.ft1 -parm NDSIZE -dim CUR_XDIM`)
set carOld = (`getParm -in test.ft1 -parm NDCAR  -dim CUR_XDIM`)
set rfp    = (`procparInfo.tcl -parm rfp -noverb`)
set rfl    = (`procparInfo.tcl -parm rfl -noverb`)
#
set refLoc = (`MATH "1 + ($size - 1)*(1 - div($rfl,$sw))"`)
#
ref1D.tcl -in test.ft1 -ref $rfp -x1 $refLoc -xn $refLoc
ref1D.tcl -in test2.ft1 -ref $rfp -x1 $refLoc -xn $refLoc
# Use the new ppm calibration to back-adjust values for the time-domain data.
set carNew = (`getParm -in test.ft1 -parm NDCAR -dim CUR_XDIM`)
echo $carOld $carNew
showhdr -verb test.fid | fgrep PPM
echo ""
ref1D.tcl -in test.fid -xAdj `MATH "$carNew - $carOld"`
showhdr -verb test.fid | fgrep PPM
echo ""
cd script
