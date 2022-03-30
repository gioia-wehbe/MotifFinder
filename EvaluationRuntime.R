#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# RuntimeEvaluation: 
# A program that gathers the p-values of SNPs belonging to motifs from vriety of datasets generated using all possible combinations of 
# parameters. Then the program summarizes the metrics of performance (accuracy, true positive, true negatives...) for all
# the datasets.
# Author: Gioia Wehbe
# LAST UPDATED: 22-02-15
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


print("Inside RuntimeEvaluation.R")
ptm <- proc.time()

print("Generating evaluation files...")

OutputText=getOutputText()

variable_parameter=getVariableParameter()     #name of the parameter that will be changing

######################################################################################################
# Performing the evaluation of the runtime
######################################################################################################

data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingDir()
}else if(data_type=="Testing"){
  data_type_dir=getTestingDir()
}else {
  data_type_dir=getCompleteDataRuntimeEvaluationDir()
}

runtime_files_dir=""
if(OutputRuns)
{
  runtime_files_dir=getRunTimeEvaluationPerFoldDir()
  if(data_type=="complete")
  {
    runtime_files_dir=getRunTimeEvaluationDir()
  }
} else{
  runtime_files_dir=getRunTimeEvaluationDir()
}
runtime_summary_file=paste(runtime_files_dir,"/Algorithm_Runtime_Summary.csv",sep="")



col_names=getOutputFileNamesShort()  
print("col_names")
print(col_names)
if(OutputPlots)
{
  col_names=getColNames()   #should return a list of values for the current variable parameter
  col_names=paste(variable_parameter,col_names,sep="_")
}
col_names=c("Code Section",col_names)

runtime_files=list.files(runtime_files_dir,full.names=TRUE,pattern="Algorithm*")

current_runtime_file_name=runtime_files[1]
runtime_table=read.csv(current_runtime_file_name,header=TRUE,sep=",",colClasses="character")

runtime_summary_table=matrix(runtime_table[,1],nrow=length(runtime_table[,1]))

for(i in 1:length(runtime_files))
{
  current_runtime_file_name=runtime_files[i]
  runtime_table=read.csv(current_runtime_file_name,header=TRUE,sep=",",colClasses="character")
  time_column=runtime_table[,2]
  runtime_summary_table=cbind(runtime_summary_table,time_column)
}

print("col_names")
print(col_names)
colnames(runtime_summary_table)=col_names

if(OutputRuns)
{
  write.csv(runtime_summary_table,runtime_summary_file,col.names=TRUE)
}

all_runs_runtime_summary_file=paste(data_type_dir,"/Table_runtime_summary_all_runs.csv",sep="")

if(file.exists(all_runs_runtime_summary_file)==FALSE) {
  write.table(runtime_summary_table,file=all_runs_runtime_summary_file,row.names=FALSE,append=TRUE,sep=",")
} else if(file.exists(all_runs_runtime_summary_file)==TRUE) {
  write.table(runtime_summary_table,file=all_runs_runtime_summary_file,row.names=FALSE,append=TRUE,sep=",",col.names=FALSE)
}

print("RuntimeEvaluation CODE ENDED!")   #print out in terminal an indication that the code has ended 

time=proc.time() - ptm  
print(paste("TIME in sec: ",time['elapsed']))
