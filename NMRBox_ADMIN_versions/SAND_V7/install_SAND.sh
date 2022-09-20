#install SAND V7 08-30-2022 OS
#This software is provided "as is" with no warranties of any kind,  and without liability for use or loss of use.
#need to define environment variable "SORIGIN" or in each .sh files: refer to the curret directory with libraries. Also PATH should include this current directory as well.

svn export https://github.com/edisonomics/SAND/trunk/scripts/NMRBox/pipe_scripts
svn export https://github.com/edisonomics/SAND/trunk/src
git clone https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA.git
#mkdir data temp log submissions data_preprocessed
chmod a+rx ./PRESAND.sh
chmod a+rx ./SAND.sh
chmod a+rx ./SUBMITSAND.sh
chmod a+rx ./pipe_scripts/*

#export PATH=$PATH:/home/nmrbox/osanati/read-only/test1/

#export SANDPATH=/home/nmrbox/osanati/read-only/test1
#export PRESAND=$SANDPATH/preprocess.sh
#export SAND=$SANDPATH/sand.sh
#export SUBMITSAND=$SANDPATH/submitsand.sh

#alias PRESAND=$SANDPATH/preprocess.sh
#alias SAND=$SANDPATH/sand.sh
#alias SUBMITSAND=$SANDPATH/submitsand.sh 
