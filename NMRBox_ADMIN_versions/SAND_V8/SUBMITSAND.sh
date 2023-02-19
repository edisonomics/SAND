#!/bin/bash

# submitsand V8 2/5/2023 frank.delaglio@nist.gov
# submitsand V7 12-28-2022 OS MODIFIED to include LAMDAMIN, SEED_12/28/2022 OS
#
#  This software is provided "as is" with no warranties of any kind,
#  and without liability for use or loss of use.

#SANDPATH:  where libraries are
#USERF:     where user put experiment folders
#USERL:     where last_preprocessed folder = WDIR
#WDIR:      same as USERL
#ORIGIN:    where commands are called
#MATLAB_EXE MATLAB command

SUBMITSAND_VERSION=V8-02_05_2023
PREF_MATLAB_EXE=/usr/software/MATLAB/R2022a/bin/matlab
MATLAB_EXE=Auto
ORIGIN=$(pwd)

echo ""
echo "SAND Spectral Decomposition"
echo "Yue Wu, Omid Sanati, Mario Uchimiya, Krish Krishnamurthy, Art Edison, and Frank Delaglio"
echo ""
echo "SUBMITSAND.sh SUBMITSAND Version:    ${SUBMITSAND_VERSION}"
echo "SUBMITSAND.sh Libraries and Scripts: $SANDPATH"
echo ""

export USERF=$ORIGIN
export NSIG=7    #max number of peaks in each bin
export NTRY=2000 #number of iteration 
export TEMP=80   #temperature in MCMC
export LAMBDAMAX=15
export LAMBDAMIN=0.01
export SEED=1

in1=0
in2=0
cores=40
name=0
MASH=0
h1=0
okinput=0
okinput2=0

while getopts a:b:x:O:z:s:n:t:I:L:K:D:M:h flag
do
    case "${flag}" in
        a) in1=${OPTARG};;	#a; low
        b) in2=${OPTARG};;	#b; high
        x) MASH=${OPTARG};;	#x; server to submit job
        O) name=${OPTARG}	#O; name of job
        	okinput2=1;;
        z) cores=${OPTARG};;	#z; number of cores
        s) NSIG=${OPTARG};;	#s; max number of peaks in each bin
        n) NTRY=${OPTARG};;	#i; number of iteration 
        t) TEMP=${OPTARG};;	#t; temparature in MCMC
        I) USERF_add=${OPTARG}
        	okinput=1
        	USERF=$ORIGIN/${OPTARG};;        
        L) LAMBDAMAX=${OPTARG};; 
        K) LAMBDAMIN=${OPTARG};;
        D) SEED=${OPTARG};;  
        M) MATLAB_EXE=${OPTARG};;
        h) h1=1;;      

    esac
done

cd $USERF
cd ..

export USERL=$ORIGIN/SAND_last_run  # write stuff there and then copy to their folders

cd $USERL

MATLAB_EXE=$(sand1D.com -matlab $MATLAB_EXE)

count=($(ls -d ./res/nmrpipe_dir/*/))

if (( $h1 == 1 )); then
	echo "This command uses HT-CONDOR distributed computing to perform SAND on the preprocessed data"
	echo "It has following options"
	echo "  -x: Server for SAND execution."
	echo "  -O: Output directory."
	echo "  -I: Input directory."
	echo "  -a  First experiment ID to process."
	echo "  -b: Last experiment ID to process."
	echo "  -z: Number of cores for each spectrum."
	echo "  -s: Maximum number of signals allowed in each bin."
	echo "  -n: Number of iterations in the Markov chain Monte Carlo (MCMC) decomposition process."
	echo "  -t: Temperature for MCMC."
	echo "  -L: Max decay, Hz."
	echo "  -K: Min decay, Hz."
	echo "  -D: Random number seed value."
	echo "  -M: MATLAB Executable."
	echo "Default values:"
	echo "  -z 100"
	echo "  -a 1 -b (total number of preprocessed folders)"
	echo "  -s 7"
	echo "  -n 2000"
	echo "  -t 80"
	echo "  -L 15"
	echo "  -K 0.01"
	echo "  -D 1"
	echo "  -M ${MATLAB_EXE} (Auto Mode)"
	echo "  -M ${PREF_MATLAB_EXE} (Preferred)"
	exit 1
fi

echo "SUBMITSAND.sh Working Directory:     $ORIGIN"
echo "SUBMITSAND.sh Last Run Folder:       $USERL"
echo "SUBMITSAND.sh User Output:           $USERF"
echo "SUBMITSAND.sh MATLAB Version:        $MATLAB_EXE"

if (( $okinput != 1 )); then
	echo "SUBMITSAND.sh Argument Error: specify the input folder with the -I option."
	exit 1
fi
if (( $okinput2 != 1 )); then
	echo "SUBMITSAND.sh Argument Error: specify the output folder with the -O option."
	exit 1
fi

t=${#count[*]}
((t += -1)) #length of all preprocessed data under ./res/nmrpipe_dir/
#echo t=$t
if (( $in1 == 0 && $in2 == 0 )); then
in1=1
in2=$t
elif (( $in1 > $t || $in2 > $t || $in1 < 1 )); then
	echo $in1 $in2 $t 
	echo "SUBMITSAND.sh Argument Error: bad experiment ID."
	exit 1
fi

if [[ ${MASH} == 0 || ${name} == 0 ]]; then
	echo "SUBMITSAND.sh Argument Error: bad server name."
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

allnames=$(ls ./log)
for i in ${allnames}
 do
	if [[ ${name} == $i ]]; then
	        echo "SUBMITSAND.sh Argument Error: Job name ${name} already exists."
		exit 1
	fi
 done

#echo ${MASH}
allinput='all'

if [[ ${MASH} != ${allinput} ]]; then

#list of available servers
rm -rf log/${name}
mkdir ./log/${name}
mkdir ./log/${name}/condor/
mkdir ./log/${name}/results 
#MASH='babbage.nmrbox.org'

#create submission files
cd ..
cat << EOF > ${USERL}/submissions/${name}.sub
universe = vanilla
executable = $SANDPATH/SAND.sh 
arguments = -a ${in1} -b ${in2} -s ${NSIG} -n ${NTRY} -t ${TEMP} -O ${name} -z ${cores} -c 1 -I ${USERF_add}  -L ${LAMBDAMAX} -K ${LAMBDAMIN} -D ${SEED} -M ${MATLAB_EXE}
request_memory = 5G
request_cpus = ${cores}
request_gpus = 0
getenv=True
requirements = (Machine==${MASH})
#requirements = (Target.Production == True)
#+Production = True
output = ${USERL}/log/${name}/condor/${name}.stdout
error = ${USERL}/log/${name}/condor/${name}.stderr
log = ${USERL}/log/${name}/condor/${name}.log
notification = Error
# Rely on a shared filesystem
transfer_executable = FALSE
should_transfer_files = NO
#JobName = ${name}
queue
EOF


#submit the job

echo "SUBMITSAND.sh CONDOR Submission:     condor_submit ${USERL}/submissions/${name}.sub"
condor_submit ${USERL}/submissions/${name}.sub

echo "SUBMITSAND.sh Job Name:          ${name}" 
echo "SUNMITSAND.sh Experiment ID:     $in1 to $in2."
echo "SUBMITSAND.sh Server:            ${MASH} with ${cores} cores." 
echo ""
echo "SUBMITSAND.sh Check the job status with 'condor_q'"

fi

if [[ ${MASH} == ${allinput} ]]; then    ####### ALL MACHINES

name2=1
name2=(${name}${name2})
allnames=$(ls ./log)
for i in ${allnames}
 do
	if [[ ${name2} == $i ]]; then
	        echo "SUBMITSAND.sh Argument Error: Job name ${name2} already exists."
		exit 1
	fi
 done

(( in1_2 = ${in1} ))
(( len = $in2 - $in1 ))
(( len += 1 ))
i=0
for i in $(seq 1 1 $len)   
do
cd $USERL
#list of available servers

name_a=(${name}${i})
rm -rf log/${name_a}
mkdir ./log/${name_a}
mkdir ./log/${name_a}/condor/
mkdir ./log/${name_a}/results 
#create submission files
cd ..
cat << EOF > ${USERL}/submissions/${name_a}.sub
universe = vanilla
executable = $SANDPATH/SAND.sh 
arguments = -a ${in1_2} -b ${in1_2} -s ${NSIG} -n ${NTRY} -t ${TEMP} -O ${name_a} -z ${cores} -c 1 -I ${USERF_add} -L ${LAMBDAMAX} -K ${LAMBDAMIN} -D ${SEED} -M ${MATLAB_EXE}
request_memory = 5G
request_cpus = ${cores}
request_gpus = 0
getenv=True
#
# fd 2/11/2023: requirements should be SAND version >= 8, NMRPipe Version >= 11.4
requirements = ((Target.Release == "2022.33") ||  (Target.Release == "2022.41") ||  (Target.Release == "2022.43") ||  (Target.Release == "2022.45"))
#requirements = (Target.Production == True)
#ASK JON about this:
#+Production = True
output = ${USERL}/log/${name_a}/condor/${name_a}.stdout
error = ${USERL}/log/${name_a}/condor/${name_a}.stderr
log = ${USERL}/log/${name_a}/condor/${name_a}.log
notification = Error
# Rely on a shared filesystem
#ASK JON about this:
transfer_executable = FALSE
should_transfer_files = NO
#JobName = ${name_a}
queue
EOF

#submit the job

echo "SUBMITSAND.sh CONDOR Submission: condor_submit ${USERL}/submissions/${name_a}.sub"
condor_submit ${USERL}/submissions/${name_a}.sub

echo "SUBMITSAND.sh Job Name:          ${name_a}" 
echo "SUNMITSAND.sh Experiment ID:     $in1_2."
echo "SUBMITSAND.sh Server:            Submitted to pool with ${cores} cores." 
echo ""
echo "SUBMITSAND.sh Check the job status with 'condor_q'"

(( in1_2 += 1 ))

done
fi

