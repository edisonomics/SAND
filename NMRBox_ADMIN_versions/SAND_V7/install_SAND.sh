#install SAND V7 09-20-2022 OS
#!/bin/bash
ORIGIN=$(dirname $(readlink -f $0))
#This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.

svn export https://github.com/edisonomics/SAND/trunk/scripts/NMRBox/pipe_scripts
svn export https://github.com/edisonomics/SAND/trunk/src
git clone --depth 1 https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA.git
#mkdir data temp log submissions data_preprocessed
cp $ORIGIN/PRESAND.sh .
cp $ORIGIN/SAND.sh .
#chmod a+rx ./SUBMITSAND.sh
chmod a+rx ./pipe_scripts/*

#export PATH=$PATH:/home/nmrbox/osanati/read-only/test1/
#export SANDPATH=/home/nmrbox/osanati/read-only/test1

