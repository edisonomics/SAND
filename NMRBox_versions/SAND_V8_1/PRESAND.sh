#!/bin/bash
# updated 2/13/2023 fd
# preprocess SAND V8: 1/11/2023 frank.delaglio@nist.gov
# This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.

#export SANDPATH=/public/groups/sand/fdelaglio/sandV8

PRESAND_VERSION=V8-02_06_2023
PREF_MATLAB_EXE=/usr/software/MATLAB/R2022a/bin/matlab
MATLAB_EXE=Auto
ORIGIN=$(pwd)

echo ""
echo "SAND Spectral Decomposition."
echo "Yue Wu, Omid Sanati, Mario Uchimiya, Krish Krishnamurthy, Art Edison, and Frank Delaglio"
echo ""
echo "PRESAND.sh PRESAND.sh Version:    ${PRESAND_VERSION}"
echo "PRESAND.sh Libraries and Scripts: $SANDPATH"
echo ""

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

while getopts p:l:u:g:M:I:h flag
do
    case "${flag}" in
        p) PREP=${OPTARG};;	  #p; already preprocessed(1) or not(0)
        l) LPPM=${OPTARG};;	  #l; lowest ppm
        u) UPPM=${OPTARG};;	  #c; highest ppm
        g) GBIN=${OPTARG};;	  #g; number of bins to be grouped
	M) MATLAB_EXE=${OPTARG};; #M; path to matlab
        I) USERF=$ORIGIN/${OPTARG}
        	okinput=1;;
        h) h1=1;;

     
    esac
done

MATLAB_EXE=$(sand1D.com -matlab $MATLAB_EXE)

export WDIR=$USERF

if (( $h1 == 1 )); then
	echo "PRESAND.sh performs preprocessing on the data to prepare it for SAND:"
	echo "  -I Input folder with folders 1 2 3 ... with test.ft in each."
	echo "  -u Upper PPM limit of analysis region."
	echo "  -l Lower PPM limit of analysis region."
	echo "  -g Number of bins for each each cluster."
	echo "  -p 1=Data is already preprocessed, 0=Data is Not Preprocessed."
	echo "  -M Path to MATLAB Command (Not all matlab versions are supported)."
	echo "Default values:"
	echo "  -p 1 (Data is already preprocessed)"
        echo "  -u 10.0"
	echo "  -l -0.4"
	echo "  -g 3"
	echo "  -M ${MATLAB_EXE} (Auto Mode)"
	echo "  -M ${PREF_MATLAB_EXE} (Preferred)"
	exit 1
fi

if (( $okinput != 1 )); then
	echo "PRESAND.sh Argument Error. Specify a valid input directory with -I option."
	exit 1
fi

if (( $PREP != 0 && $PREP != 1 )); then
	echo "PRESAND.sh Argument Error. Option -p (preprocessed status) should be 1 or 0."
	exit 1
fi

if (( $(echo "$LPPM > $UPPM" |bc -l) )); then
	fTmp=$LPPM
	LPPM=$UPPM
	UPPM=$fTmp
fi

if [[ $GBIN < 1 ]]; then
	echo "PRESAND.sh Argument Error. Number of bins must be greater than 0."
	exit 1
fi

#
# fd 2/13/2023 to minimize un-needed user support tickets about the notice:

# if [[ "$MATLAB_EXE" != "$PREF_MATLAB_EXE" ]]; then
#    echo ""
#    echo "SAND.sh Notice. Using a MATLAB Version other than the preferred one."
#    echo "SAND.sh Notice. In Use:    $MATLAB_EXE"
#    echo "SAND.sh Notice. Preferred: $PREF_MATLAB_EXE"
#    echo ""
# fi

export USERL=$ORIGIN/SAND_last_run

rm -f -r ${USERL}
mkdir ${USERL} ${USERL}/log ${USERL}/submissions ${USERL}/temp
cd $USERL
cd ..

echo "PRESAND.sh Working Directory:     $ORIGIN"
echo "PRESAND.sh Last Run Folder:       $USERL"
echo "PRESAND.sh User Output:           $USERF"
echo "PRESAND.sh MATLAB Version:        $MATLAB_EXE"

cd $USERF
cd ..
export USERF=$USERF

echo ""
echo "PRESAND.sh Preprocessing (PID=$$) starting."
echo "PRESAND.sh Several MATLAB windows will appear and close."
time ${MATLAB_EXE} < $SANDPATH/PRESAND.m > $USERL/log/preprocess$$.log 2> ${USERL}/log/preprocess$$.err
# -nojvm -nosplash -nodisplay
