#!/bin/bash
# preprocess SAND V7: 09-19-2022 OS
#This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.


#export SANDPATH=/home/nmrbox/osanati/read-only/test1
echo "where libraries are: $SANDPATH"
ORIGIN=$(pwd)
echo "where it is called: $ORIGIN"
export TDD=$SANDPATH/src 
export ELSM=$SANDPATH/Edison_Lab_Shared_Metabolomics_UGA
export NCPUS=$(nproc)
export PREP=1
export LPPM=-0.4
export UPPM=10 
export GBIN=3
export USERF=$ORIGIN
h1=0
okinput=0
while getopts p:l:u:g:I:h flag
do
    case "${flag}" in
        p) PREP=${OPTARG};;	#p; already preprocessed(1) or not(0)
        l) LPPM=${OPTARG};;	#l; lowest ppm
        u) UPPM=${OPTARG};;	#c; highest ppm
        g) GBIN=${OPTARG};;	#g; number of bins to be grouped
        I) USERF=$ORIGIN/${OPTARG}
        	okinput=1;;
        h) h1=1;;

     
    esac
done

export WDIR=$USERF

if (( $h1 == 1 )); then
	echo "this command performs preprocessing on the data to prepare it for SAND"
	echo "it has following options"
	echo "	-I input folder with bruker folders inside with test.ft in each"
	echo "  -l and -u: lower/ upper PPM range for preprocessed data"
	echo "  -g: number of bins in each cluster"
	echo "default values:"
	echo "  -p 1 (data is already preprocessed)"
	echo "  -l -0.4 u 10"
	echo "  -g 3"
	exit 1
fi

if (( $okinput != 1 )); then
	echo "try again! wrong inputs: choose the input folder with -I!"
	exit 1
fi

if (( $PREP != 0 && $PREP != 1 )); then
	echo "try again! wrong inputs: choose you preprocessed status by 1 or 0!"
	exit 1
fi

if [[ $LPPM > $UPPM ]]; then
	echo "try again! wrong inputs: choose lower/upper limit ppms correctly!"
	exit 1
fi
if [[ $GBIN < 1 ]]; then
	echo "try again! wrong inputs: choose number of bins in each group correctly!"
	exit 1
fi

export USERL=$ORIGIN/SAND_last_run  # write stuff there and then copy to their folders
rm -f -r ${USERL}
mkdir ${USERL} ${USERL}/log ${USERL}/submissions ${USERL}/temp
cd $USERL
cd ..
echo "where is last run folder: $USERL"

echo "where user refers: $USERF" # data should be there
cd $USERF
cd ..
export USERF=$USERF

echo "Preprocessing (PID=$$) started"
time /usr/software/MATLAB/R2022a/bin/matlab < $SANDPATH/PRESAND.m > $USERL/log/preprocess$$.log 2> ${USERL}/log/preprocess$$.err
#-nojvm -nosplash -nodisplay  > tda$$.log 2>tda$$.err
