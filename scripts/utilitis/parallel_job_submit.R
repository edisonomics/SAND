# batch submision script for spectral decomposition
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
dir=""
shelltempt=paste0(dir,"submit.sh")
matlabtempt=paste0(dir,"test_region_separa_deconv_hpc.m")
tabrun=read.table("runtab.txt",header=TRUE)
for(runi in 1:dim(tabrun)[1]){
  runmatlab=paste0("test_region_separa_deconv_hpc",runi,".m")
  runshell=paste0("runshell",runi,".sh")
  # change shell files
  system(paste0("cp ",shelltempt," ",runshell))
  lines=readLines(runshell)
  chline_ind=str_which(string=lines,pattern="^time")
  lines[chline_ind]=paste0("time matlab -nodisplay < ",runmatlab)
  cat(lines,file=runshell,sep="\n")
  # change matlab files
  system(paste0("cp ",matlabtempt," ",runmatlab))
  lines=readLines(runmatlab)
  for(col in colnames(tabrun)){
    chline_ind=str_which(string=lines,pattern=paste0("^",col,"_arg="))
    val=tabrun[runi,col]
    if(col=="dataset"){
      val=paste0("'",val,"'");
    }
    lines[chline_ind]=paste0(col,"_arg=",val)
  }
  cat(lines,file=runmatlab,sep="\n")
  #
  submitcommand=paste0("sbatch ",runshell)
  print(submitcommand)
  system(submitcommand)
}
