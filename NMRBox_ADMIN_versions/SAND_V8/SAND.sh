#!/bin/bash
# SAND V8: 2/13/2023 frank.delaglio@nist.gov
# This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.
# MODIFIED to include LAMDAMIN, SEED, NCPU_12/28/2022 OS

#SANDPATH:  where libraries are
#USERF:     where user put experiment folders
#USERL:     where last_preprocessed folder = WDIR
#WDIR:      same as USERL
#ORIGIN:    where commands are called
#MATLAB_EXE MATLAB command

SAND_VERSION=V8-02_06_2023
PREF_MATLAB_EXE=/usr/software/MATLAB/R2022a/bin/matlab
MATLAB_EXE=Auto
ORIGIN=$(pwd)

#export SANDPATH=/public/groups/sand/fdelaglio/sandV8

echo ""
echo "SAND Spectral Decomposition"
echo "Yue Wu, Omid Sanati, Mario Uchimiya, Krish Krishnamurthy, Art Edison, and Frank Delaglio"
echo ""
echo "SAND.sh SAND.sh Version:       ${SAND_VERSION}"
echo "SAND.sh Libraries and Scripts: $SANDPATH"
echo ""

export USERF=$ORIGIN
export TDD=$SANDPATH/src 
export ELSM=$SANDPATH/Edison_Lab_Shared_Metabolomics_UGA
export NCPUS=0
export NSAMPLE=1
export NSIG=7    #max number of peaks in each bin
export NTRY=2000 #number of iteration 
export TEMP=80   #temperature in MCMC
export LAMBDAMAX=15
export LAMBDAMIN=0.01
export SEED=1

in1=0
in2=0
cluster=0
name=0
h1=0
okinput1=0
okinput2=0

while getopts a:b:c:O:s:n:t:z:I:L:K:D:C:M:h flag
do
    case "${flag}" in
        a) in1=${OPTARG};;	#a; low
        b) in2=${OPTARG};;	#b; high
        c) cluster=${OPTARG};;	#c; cluster or not
        O) name=${OPTARG}
        	okinput2=1;;	#O; name of job
        s) NSIG=${OPTARG};;
        n) NTRY=${OPTARG};;
        t) TEMP=${OPTARG};;
        z) NCPUS=${OPTARG};;
        I) USERF=$ORIGIN/${OPTARG}
        	okinput1=1;;
        L) LAMBDAMAX=${OPTARG};;
        K) LAMBDAMIN=${OPTARG};;
        D) SEED=${OPTARG};;
        C) NCPUS=${OPTARG};;
        M) MATLAB_EXE=${OPTARG};;
        h) h1=1;;    

    esac
done

MATLAB_EXE=$(sand1D.com -matlab $MATLAB_EXE)
maxCores=$(nproc)

if [[ $NCPUS == 0 ]]; then
   NCPUS=$maxCores
#   NCPUS=$(IMATH "1 + ($NCPUS/2)")
	NCPUS=$((NCPUS/2))
fi
   
if (( $NCPUS > $maxCores )); then
   NCPUS=$maxCores
fi

cd $USERF
cd ..

export USERL=$ORIGIN/SAND_last_run  # write stuff there and then copy to their folders
export WDIR=$USERL
export TEMPDIR=$USERL/temp
export LOGDIR=$USERL/log

if (( $h1 == 1 )); then
	echo "Perform SAND decomposition on data preprocessed with PRESAND:"
	echo "  -O: Output directory."
	echo "  -a: First experiment ID to process."
	echo "  -b: Last experiment ID to process."
	echo "  -s: Max number of signals allowed in each bin."
	echo "  -n: Number of iterations in the Markov chain Monte Carlo (MCMC) decomposition process."
	echo "  -t: Temperature for MCMC."
	echo "  -I: Input directory."
	echo "  -L: Max linewidth, Hz."
	echo "  -K: Min linewidth, Hz."
	echo "  -D: Random number seed value."
	echo "  -C: Number of cores for interactive mode."
	echo "Default values:"
	echo "  -a 1 -b (total number of preprocessed spectra)"
	echo "  -s 7"
	echo "  -n 2000"
	echo "  -t 80"
	echo "  -L 15"
	echo "  -K 0.01"
	echo "  -D 1"
	echo "  -C (max_cores)/2"
	echo "  -M ${MATLAB_EXE} (Auto Mode)"
	echo "  -M ${PREF_MATLAB_EXE} (Preferred)"
	
	exit 1
fi

echo "SAND.sh Working Directory:     $ORIGIN"
echo "SAND.sh Last Run Folder:       $USERL"
echo "SAND.sh User Output:           $USERF"
echo "SAND.sh MATLAB Version:        $MATLAB_EXE"

if (( $okinput1 != 1 )); then
	echo "SAND.sh Argument Error. Specify a valid input directory with the -I option."
	exit 1
fi
if (( $okinput2 != 1 )); then
	echo "SAND.sh Argument Error. Specify an output directory with the -O option."
	exit 1
fi

#
# fd 2/13/2023 to minimize un-needed user support tickets about the notice:

# if [[ "$MATLAB_EXE" != "$PREF_MATLAB_EXE" ]]; then
#   echo ""
#   echo "SAND.sh Notice. Using a MATLAB Version other than the preferred one."
#   echo "SAND.sh Notice. In Use:    $MATLAB_EXE"
#   echo "SAND.sh Notice. Preferred: $PREF_MATLAB_EXE"
#   echo ""
#fi

cd $USERL

count=($(ls -d ./res/nmrpipe_dir/*/))
export JOBNAME=${name}

t=${#count[*]}
((t += -1)) #length of all preprocessed data under ./res/nmrpipe_dir/
#echo $t
if (( $in1 == 0 && $in2 == 0 )); then
in1=1
in2=$t
elif (( $in1 > $t || $in2 > $t || $in1 < 1 )); then
	echo "SAND.sh Argument Error. An incorrect experiment number was specified."
	exit 1
fi

if [[ ${cluster} == 0 ]]; then
	if [[ ${name} == 0 ]]; then
		echo "SAND.sh Argument Error. Error finding input ${name}."
		exit 1
	fi

	allnames=$(ls ./log)
	for i in ${allnames}
	 do
		if [[ ${name} == $i ]]; then
			echo "SAND.sh Argument Error. Job name ${name} already exists."
			exit 1
		fi
	 done
	 
        rm -rf log/${name}
	mkdir ./log/${name}
	mkdir ./log/${name}/results 
	 
fi

NSAMPLE=$in1   
THISDIR=$(pwd)

echo ""
echo "SAND.sh Time-domain analysis (process name= ${name}) starting on $NCPUS threads."
echo "SAND.sh Experiment ID: $in1 to $in2."

while [ $NSAMPLE -ge $in1 ] && [ $NSAMPLE -le $in2 ]
do
	echo ""
	echo "SAND.sh Analyzing Experiment ID $NSAMPLE"
	time ${MATLAB_EXE} < $SANDPATH/SAND.m > $LOGDIR/${name}/process${name}.log 2> $LOGDIR/${name}/process${name}.err
	echo "SAND.sh Finished Experiment ID $NSAMPLE"

	if [[ -f "$ORIGIN/pipeUpdate.com" ]]; then
	   echo "SAND.sh Updating NMRPipe Files for Experiment ID $NSAMPLE ..."
	   echo ""
	   cd $ORIGIN
	      pipeUpdate.com -firstID $NSAMPLE -lastID $NSAMPLE
	   cd $THISDIR
	   echo ""
	fi

	((NSAMPLE += 1))
done

