####################################################
# Runs SPLSDA for Benchmarking
# last updated: 22-06-15
####################################################

print("inside RunSPLSDAforTiming.R")

current_fold = getCurrentFold()
if(data_type=="complete")
{
  print("inside if(data_type==complete)")
  current_fold=getData() #evaluate time on complete data set
  print("dim(current_fold)")
  print(dim(current_fold))
  
}
parameter_table = read.table(paste(getwd(),"/Parameters/parameters.txt",sep=""),header=FALSE,sep="\t",strip.white=TRUE)
parameter_table=as.matrix(parameter_table)  #import input file containing parameters as a table and convert it to a matrix
num_classes=as.numeric(parameter_table[1,2]) 
print("num_classes")
print(num_classes)
num_SNPs= as.numeric(parameter_table[3,2]) 

output_file_name=parameter_table[7,2]   #output files name base on input parameters

output_file_name_list=unlist(strsplit(output_file_name,"_"))
print("output_file_name_list")
print(output_file_name_list)
MperC_index=length(output_file_name_list)-4
print("MperC_index")
print(MperC_index)
str_MperC=output_file_name_list[MperC_index]
print("str_MperC")
print(str_MperC)
LM_index=length(output_file_name_list)-3
print("LM_index")
print(LM_index)
str_LM=output_file_name_list[LM_index]
print("str_LM")
print(str_LM)

num_motifs_per_class=as.numeric(unlist(regmatches(str_MperC,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",str_MperC))))
print("num_motifs_per_class")
print(num_motifs_per_class)
length_motif=as.numeric(unlist(regmatches(str_LM,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",str_LM))))
print("length_motif")
print(length_motif)

colnames(current_fold)=append(1:num_SNPs,"Class")


class_labels=current_fold[,'Class']
X=current_fold[,1:num_SNPs]
Y=as.factor(class_labels)

ncomp_var=num_motifs_per_class*num_classes               #Gioia: decide on the best value for the number of components (motifs_per_class x num_classes)
print("ncomp_var")
print(ncomp_var)
keepX_var=rep(length_motif,ncomp_var)        #Gioia: decide on the best value for the number of variables to keep for each component (Length_of_motif)
print("keepX_var")
print(keepX_var)

timing_summary_table_CPU=matrix(ncol=2)
colnames(timing_summary_table_CPU)=c("Code Section","Time in Seconds")


start=proc.time()
data_splsda <- splsda(X, Y, ncomp = ncomp_var, keepX = keepX_var)

splsda_duration=proc.time()-start

splsda_duration_CPU=splsda_duration['user.self']+splsda_duration['sys.self']
splsda_duration_CPU=round(splsda_duration_CPU,digits=3)
splsda_duration_CPU_row=c("Duration for running SPLSDA Feature Selection Algorithm",splsda_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,splsda_duration_CPU_row)

print("splsda_duration")
print(splsda_duration)


start_select_var=proc.time()
all_selected_vars_names=c()

for(i in 1:ncomp_var)
{
  current_comp=i
  current_comp_selected_vars=selectVar(data_splsda,comp=current_comp)
  current_comp_selected_vars_names=current_comp_selected_vars$name
  all_selected_vars_names=c(all_selected_vars_names,current_comp_selected_vars_names)
}

all_selected_vars_names=unique(all_selected_vars_names)

duration_select_var=proc.time()-start_select_var

duration_select_var_CPU=duration_select_var['user.self']+duration_select_var['sys.self']
duration_select_var_CPU=round(duration_select_var_CPU,digits=3)
duration_select_var_CPU_row=c("Duration for selecting Vars",duration_select_var_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,duration_select_var_CPU_row)

print("duration_select_var")
print(duration_select_var)

total_splsda_duration=duration_select_var+splsda_duration

total_splsda_duration_CPU=total_splsda_duration['user.self']+total_splsda_duration['sys.self']
total_splsda_duration_CPU=round(total_splsda_duration_CPU,digits=3)
total_splsda_duration_CPU_row=c("Total Duration for SPLSDA Feature Selection Algorithm",total_splsda_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,total_splsda_duration_CPU_row)

timing_summary_table_CPU=timing_summary_table_CPU[-1,]



# all_selected_vars_names=sort(as.numeric(all_selected_vars_names))
# all_snps=colnames(X)
# all_selected_vars_boolean=all_snps %in% all_selected_vars_names

# print("current_output_file_name_short")
# print(current_output_file_name_short)

# all_selected_vars_names_output_file=paste(getSplsdaPredictedDiffLocations(),"/SPLSDA_Diff_locations_",current_output_file_name_short
#                                           ,".txt",sep="")
# all_selected_vars_boolean_output_file=paste(getSplsdaPredictedDiffLocations(),"/SPLSDA_Diff_locations_boolean_",current_output_file_name_short
#                                             ,".txt",sep="")
# 
# write(all_selected_vars_names,all_selected_vars_names_output_file,ncolumns=1)
# write(all_selected_vars_boolean,all_selected_vars_boolean_output_file,ncolumns=1)

# print("all_selected_vars_names_output_file")
# print(all_selected_vars_names_output_file)
# 
# print("all_selected_vars_names")
# print(all_selected_vars_names)

time_file_name=""
if(OutputRuns)
{
  time_file_name=paste(getRunTimeEvaluationPerFoldDir(),"/SplsdaRuntime_",current_output_file_name_short,".csv",sep="")
  
  if(data_type=="complete")
  {
    time_file_name=paste(getRunTimeEvaluationDir(),"/SplsdaRuntime_",current_output_file_name_short,".csv",sep="")
  }
}else{
  time_file_name=paste(getRunTimeEvaluationDir(),"/SplsdaRuntime_",current_output_file_name_short,".csv",sep="")
}
write.csv(timing_summary_table_CPU,time_file_name,row.names=FALSE)

print("RunSPLSDAforTiming.R code ended")




