#!/bin/bash

#submitsand V7 09-19-2022 OS
#This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.

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

export NSIG=7  #max number of peaks in each bin
export NTRY=2000 #number of iteration 
export TEMP=80 #temperature in MCMC

in1=0
in2=0
cores=100
name=0
MASH=0
h1=0
okinput=0
okinput2=0
while getopts a:b:x:O:z:s:n:t:I:h flag
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
        h) h1=1;;

    esac
done

echo "where user refers: $USERF"
cd $USERF
cd ..

export USERL=$ORIGIN/SAND_last_run  # write stuff there and then copy to their folders
echo "where last preprocessed data is: $USERL"
cd $USERL
count=($(ls -d ./res/nmrpipe_dir/*/))

if (( $h1 == 1 )); then
	echo "this command performs SAND on the preprocessed data"
	echo "it has following options"
	echo "  -x: you must define which server to run SAND on"
	echo "  -O: you must assign a unique name to each run of SAND"
	echo "  -a and -b: you can define which experiments (adjusted numbers) to run SAND on them"
	echo "  -z: number of cores to run each experiment"
	echo "  -s: maximum number of signals (after the deconvolution process) in each bin"
	echo "  -n: number of iterations in the MCMC (deconvolution process)"
	echo "  -t: temperature in  the MCMC (deconvolution process)"
	echo "  -I: input folder"
	echo "default values:"
	echo "  -z 100"
	echo "  -a 1 -b (total number of preprocessed folders)"
	echo "  -s 7"
	echo "  -n 2000"
	echo "  -t 80"
	exit 1
fi

if (( $okinput != 1 )); then
	echo "try again! wrong inputs: choose the input folder with -I!"
	exit 1
fi
if (( $okinput2 != 1 )); then
	echo "try again! wrong inputs: choose the output folder with -O!"
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
	echo "try again! wrong experiment number!"
	exit 1
fi

if [[ ${MASH} == 0 || ${name} == 0 ]]; then
	echo "try again! wrong inputs: choose servers and names correctly!"
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
arguments = -a ${in1} -b ${in2} -s ${NSIG} -n ${NTRY} -t ${TEMP} -O ${name} -z ${cores} -c 1 -I ${USERF_add}
request_memory = 5G
request_cpus = ${cores}
request_gpus = 0
getenv=True
requirements = (Machine==${MASH})
requirements = (Target.Production == True)
+Production = True
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

condor_submit ${USERL}/submissions/${name}.sub
echo "Job '${name}' :experiment(s) $in1 to $in2."
echo "It is submitted to '${MASH}' with ${cores} cores." 
echo "Check the status with 'condor_q'"

fi

if [[ ${MASH} == ${allinput} ]]; then    ####### ALL MACHINES

name2=1
name2=(${name}${name2})
allnames=$(ls ./log)
for i in ${allnames}
 do
	if [[ ${name2} == $i ]]; then
		echo "try again! Job name already exists, choose different name correctly!"
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
arguments = -a ${in1_2} -b ${in1_2} -s ${NSIG} -n ${NTRY} -t ${TEMP} -O ${name_a} -z ${cores} -c 1 -I ${USERF_add}
request_memory = 5G
request_cpus = ${cores}
request_gpus = 0
getenv=True
requirements = ((Target.Release == "2022.37") ||  (Target.Release == "2022.35") ||  (Target.Release == "2022.33") ||  (Target.Release == "2022.31") ||  (Target.Release == "2022.29") || (Target.Release == "2022.35") || (Target.Release == "2022.33") || (Target.Release == "2022.27"))
#requirements = (Target.Production == True)
+Production = True
output = ${USERL}/log/${name_a}/condor/${name_a}.stdout
error = ${USERL}/log/${name_a}/condor/${name_a}.stderr
log = ${USERL}/log/${name_a}/condor/${name_a}.log
notification = Error
# Rely on a shared filesystem
transfer_executable = FALSE
should_transfer_files = NO
#JobName = ${name_a}
queue
EOF


#submit the job

condor_submit ${USERL}/submissions/${name_a}.sub
echo "Job '${name_a}' :experiment $in1_2"
echo "It is submitted to the pool with ${cores} cores." 
echo "Check the status with 'condor_q'"
(( in1_2 += 1 ))

done
fi

