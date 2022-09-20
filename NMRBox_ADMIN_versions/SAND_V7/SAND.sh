#!/bin/bash
#This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.
# process SAND V7: 09-13-2022 OS

#SANDPATH: where libraries are
#USERF: where user put experiment folders
#USERL: where last_preprocessed folder = WDIR
#WDIR: same as USERL
#ORIGIN: where commands are called

#export SANDPATH=/home/nmrbox/osanati/read-only/test1
echo "where libraries are: $SANDPATH"
ORIGIN=$(pwd)
echo "where it is called: $ORIGIN"
export USERF=$ORIGIN

export TDD=$SANDPATH/src 
export ELSM=$SANDPATH/Edison_Lab_Shared_Metabolomics_UGA
export NCPUS=$(nproc)
export NSAMPLE=1
export NSIG=7  #max number of peaks in each bin
export NTRY=2000 #number of iteration 
export TEMP=80 #temparature in MCMC

in1=0
in2=0
cluster=0
name=0
h1=0
okinput1=0
okinput2=0
while getopts a:b:c:O:s:n:t:z:I:h flag
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
        h) h1=1;;        

    esac
done

echo "where user refers: $USERF"
cd $USERF
cd ..

export USERL=$ORIGIN/SAND_last_run  # write stuff there and then copy to their folders
export WDIR=$USERL
export TEMPDIR=$USERL/temp
export LOGDIR=$USERL/log
cd $USERL

if (( $h1 == 1 )); then
	echo "this command performs SAND on the preprocessed data"
	echo "it has following options"
	echo "  -O: you must assign a unique name to each run of SAND"
	echo "  -a and -b: you can define which experiments (adjusted numbers) to run SAND on them"
	echo "  -z: number of cores to run each experiment"
	echo "  -s: maximum number of signals (after the deconvolution process) in each bin"
	echo "  -n: number of iterations in the MCMC (deconvolution process)"
	echo "  -t: temperature in  the MCMC (deconvolution process)"
	echo "  -I: input folder"
	echo "default values:"
	echo "  -z (total number of physical cores)"
	echo "  -a 1 -b (total number of preprocessed folders)"
	echo "  -s 7"
	echo "  -n 2000"
	echo "  -t 80"
	exit 1
fi

if (( $okinput1 != 1 )); then
	echo "try again! wrong inputs: choose the input folder with -I!"
	exit 1
fi
if (( $okinput2 != 1 )); then
	echo "try again! wrong inputs: choose the output folder with -O!"
	exit 1
fi

count=($(ls -d ./res/nmrpipe_dir/*/))
export JOBNAME=${name}

t=${#count[*]}
((t += -1)) #length of all preprocessed data under ./res/nmrpipe_dir/
#echo $t
if (( $in1 == 0 && $in2 == 0 )); then
in1=1
in2=$t
elif (( $in1 > $t || $in2 > $t || $in1 < 1 )); then
	echo "try again! wrong experiment number!"
	exit 1
fi

if [[ ${cluster} == 0 ]]; then
	if [[ ${name} == 0 ]]; then
		echo "try again! wrong inputs: choose the name correctly!"
		exit 1
	fi

	allnames=$(ls ./log)
	for i in ${allnames}
	 do
		if [[ ${name} == $i ]]; then
			echo "try again! Job name already exists, choose different name correctly!"
			exit 1
		fi
	 done
	 
	 rm -rf log/${name}
	mkdir ./log/${name}
	mkdir ./log/${name}/results 
	 
fi



NSAMPLE=$in1   


echo "time-domain analysis (process name= ${name}) started on $NCPUS threads"
echo "for experiment(s) $in1 up to $in2."

while [ $NSAMPLE -ge $in1 ] && [ $NSAMPLE -le $in2 ]
do
	echo "experiment:$NSAMPLE started"
	time /usr/software/MATLAB/R2022a/bin/matlab < $SANDPATH/SAND.m > $LOGDIR/${name}/process${name}.log 2> $LOGDIR/${name}/process${name}.err
	echo "experiment:$NSAMPLE finished"
	((NSAMPLE += 1))
done
