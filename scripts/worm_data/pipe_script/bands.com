#!/bin/csh

#
# Example of creating overlapping bands from a 1D spectrum.
# frank.delaglio@nist.gov Oct 15 2020
#
# This script prepares a series of overlapping cosine rolloff
# bands, and uses them to create a spectral series, each of which
# contains one band.
#
# The scheme makes use of a utility called "normSpec.tcl", a general purpose
# utility for normalizing a spectrum.  The utility is used by constructing a
# freeform list of operations to perform on the given spectrum for example
# "set region 12ppm 0.0ppm to 1.0". The band masks are created via operations
# that generate a rising or falling cosine-square ramp.
cd ../

set T_r=0.05
set inName          = (`getArgD $argv -in         test.ft1`)
set specPrefix      = (`getArgD $argv -specPrefix spec`)
set maskPrefix      = (`getArgD $argv -maskPrefix mask`)
set bandWidthPPM    = (`getArgD $argv -bw         $T_r`)#not for each ppm range but only for the initial spectral processing

if (`flagLoc $argv -help`) then
   echo "Create a Series of Band Selection Masks:"
   echo " -in         inName             [$inName]"
   echo " -specPrefix specPrefix         [$specPrefix]"
   echo " -maskPrefix maskPrefix         [$maskPrefix]"
   echo " -bw         bandWidthPPM       [$bandWidthPPM]"
   exit 0
endif

set outDir = $maskPrefix
set ftDir  = ${maskPrefix}_ft
set fidDir = ${maskPrefix}_fid

/bin/rm -rf $outDir $ftDir $fidDir
mkdir $outDir $ftDir $fidDir

#
# Prepare a mask which zeros the edges of the spectrum.
# This version of the spectrum is used as input for the
# other steps.
#
# From this point on, we can treat the spectrum as if it has
# no digital oversampling and no phase correction applied.

set size   = (`getParm -in $inName -parm NDSIZE -dim CUR_XDIM -fmt %.0f`)
set specX1 = (`pnt2spec $inName CUR_XDIM 1     ppm`)
set specXN = (`pnt2spec $inName CUR_XDIM $size ppm`)

set x2 = (`MATH "$specX1 - $bandWidthPPM"`)
set x3 = (`MATH "$specXN + $bandWidthPPM"`)

# check value
echo $x2
echo $x3

set edgeName    = ${maskPrefix}_edge.ft1
set edgeFTName  = ${specPrefix}_edge.ft1
set edgeFIDName = ${specPrefix}_edge.fid

/bin/rm -f $edgeName $edgeFTName $edgeFIDName

normSpec.tcl -in $inName -out $edgeName \
    -prog set      region 0%       100%     to 1.0 then \
          rampUp   region 0%       ${x2}ppm        then \
          set      region ${x2}ppm ${x3}ppm to 1.0 then \
          rampDown region ${x3}ppm 100%

addNMR -in1 $inName -in2 $edgeName -out $edgeFTName -mult

sethdr $edgeName -fdata FDDMXVAL 0 FDDMXFLAG 0 -xP0 0.0 -xP1 0.0

nmrPipe -in $edgeFTName \
| nmrPipe -fn HT \
| nmrPipe -fn FT -inv \
| nmrPipe -fn ZF -inv \
  -out $edgeFIDName -ov

sethdr $edgeName    -title mask_edges
sethdr $edgeFTName  -title spec_masked_edges
sethdr $edgeFIDName -title ift_masked_edges

#
# Create a series of overlapping masks.
# For each mask, multiply it by the original spectrum.
# Inverse process the masked result (to increased stability, window function is not removed).

set file='binrange.txt'
set rangel=(`awk '{print $1}' "$file"`)
set ranger=(`awk '{print $2}' "$file"`)
set len=${#rangel}

@ i = 1

set vFlag = "-verb"
rm record.txt
# nmrPrintf "FORMAT %%s %%6.4f %%6.4f %%6.4f %%6.4f\n" >> record.txt
# nmrPrintf "VARS NAME IN_LEFT IN_RIGHT OUT_LEFT OUT_RIGHT\n" >> record.txt
while ($i <= $len)
   set thisName = (`nmrPrintf "%s/test%03d.ft1" $outDir $i`)
   set ftName   = (`nmrPrintf "%s/test%03d.ft1" $ftDir  $i`)
   set fidName  = (`nmrPrintf "%s/test%03d.fid" $fidDir $i`)

   set x2 = $rangel[$i]
   set x3 = $ranger[$i]
   set rollOffWidthPPM=(`MATH "($x2 - $x3)/2"`)
   set x1 = (`MATH "$x2 + $rollOffWidthPPM"`)
   set x4 = (`MATH "$x3 - $rollOffWidthPPM"`)


   nmrPrintf "%s %s %s %6.4f ppm %6.4f ppm %6.4f ppm %6.4f ppm\n" $thisName $ftName $fidName $x1 $x2 $x3 $x4
   nmrPrintf "%s %6.4f  %6.4f  %6.4f  %6.4f\n" $fidName $x2 $x3 $x1 $x4 >> record.txt
   
   normSpec.tcl $vFlag -in $edgeFTName -out $thisName \
    -prog set      region 0%       100%     to 0.0 then \
          set      region ${x1}ppm ${x4}ppm to 1.0 then \
          rampUp   region ${x1}ppm ${x2}ppm        then \
          set      region ${x2}ppm ${x3}ppm to 1.0 then \
          rampDown region ${x3}ppm ${x4}ppm
   
   normSpec.tcl $vFlag -in $inName -out $ftName \
    -prog subtract min of region ${x1}ppm ${x4}ppm
  
   addNMR -in1 $ftName -in2 $thisName -out $ftName -mult
   
   nmrPipe -in $ftName \
   | nmrPipe -fn HT \
   | nmrPipe -fn FT -inv \
   | nmrPipe -fn ZF -inv \
     -out $fidName -ov

   sethdr $thisName -title `nmrPrintf mask%03d       $i`
   sethdr $ftName   -title `nmrPrintf spec_band%03d  $i`
   sethdr $fidName  -title `nmrPrintf ift_band%03d   $i`

   if ($i == 1) then
      echo ""
   endif
   
   @ i++
end

#
# In order to account for overlapping regions, in order to generate a complete
# spectrum by adding up each band, the result must be normalized (divided) by
# the sum of all the masks.

/bin/rm -f ${maskPrefix}_sum.ft1 ${maskPrefix}_sum_ft.ft1 ${maskPrefix}_sum_fid.ft1
 
sumNMR.com -in $outDir/* -out ${maskPrefix}_sum.ft1
sumNMR.com -in $ftDir/*  -out ${maskPrefix}_sum_ft.ft1

addNMR -in1 ${maskPrefix}_sum_ft.ft1 -in2 ${maskPrefix}_sum.ft1 -out ${maskPrefix}_sum_ft.ft1 -div

nmrPipe -in ${maskPrefix}_sum_ft.ft1 \
| nmrPipe -fn HT \
| nmrPipe -fn FT -inv \
| nmrPipe -fn ZF -inv \
   -out ${maskPrefix}_sum_ift.fid -ov

sethdr ${maskPrefix}_sum.ft1     -title mask_sum
sethdr ${maskPrefix}_sum_ft.ft1  -title bands_sum
sethdr ${maskPrefix}_sum_ift.fid -title ift_bands_sum

echo ""
echo "1. Input Spectrum:                            $inName"
echo "2. Mask for Edges:                            $edgeName"
echo "3. Masks for Bands                            ${maskPrefix}/test%03d.ft1"
echo "4. Input Spectrum, Edges Masked (1)x(2):      $edgeFTName"
echo "5. Inverse Processing of (4):                 $edgeFIDName"
echo "6. Spectral Bands (1)x(2)x(3):                ${maskPrefix}_ft/test%03d.fid"
echo "7. Sum of All Masks (3):                      ${maskPrefix}_sum.ft1"
echo "8. Sum of Masked Spectra (6), Divided by (7): ${maskPrefix}_sum_ft.ft1"
echo "9. Inverse Processing of (8):                 ${maskPrefix}_sum_ift.fid"

cd ./script
