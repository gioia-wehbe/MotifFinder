
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# DiffAreasEvaluation: 
# A program that performs the evaluation of DiffAreas algorithm and outputs evaluation tables as well as boxplots
# if OutputText=TRUE
# Author: Gioia Wehbe
# LAST UPDATED: 23-02-15
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print("Inside SplsdaEvaluation.R")
ptm <- proc.time()

print("Generating evaluation files...")

OutputText=getOutputText()
OutputPlots=getOutputPlots()
variable_parameter=getVariableParameter()     #name of the parameter that will be changing

#Directories of the files that will be needed 
actual_motifs_locations_dir=getActualMotifsDir()
# run_dir=getRunsDir()

training_dir=getTrainingDir()
print("training_dir")
print(training_dir)
testing_dir=getTestingDir()
print("testing_dir")
print(testing_dir)
eval_output_dir=getEvaluationOutputDir()

current_run_name=getCurrentRunName()
current_fold_name=getCurrentFoldName()

#output tables initialization, used to store the evaluation output of the differentiable areas identifier 
output_table_metrics=matrix(nrow=16)
output_table_counts=matrix(nrow=16)

#output files names where the the evaluation tables of the differentiable areas identifiers per run will be stored
splsda_differentiable_locations_rates=paste(getPredicionEvaluationDir(training_dir),"/splsda_differentiable_locations_per_allele_SEN_SPE_PRE_NPV_ACC_summary.csv",sep="")
splsda_differentiable_locations_counts=paste(getPredicionEvaluationDir(training_dir),"/splsda_differentiable_locations_per_allele_TP_FP_TN_FN_ACC_summary.csv",sep="")

#output files names used to store the evaluation of the differentiable areas identifier against all runs
evaluation_table_all_runs=paste(eval_output_dir,"/Table_splsda_differentiable_locations_per_allele_metrics_all_runs.csv",sep="")
evaluation_counts_table_all_runs=paste(eval_output_dir,"/Table_splsda_differentiable_locations_per_allele_counts_all_runs.csv",sep="")


#Input files
actual_motifs_filenames <- list.files(actual_motifs_locations_dir, full.names=TRUE)
output_file_names_short=getOutputFileNamesShort()
print("output_file_names_short")
print(output_file_names_short)
output_file_names_long=getOutputFileNamesLong()
print("output_file_names_long")
print(output_file_names_long)

######################################################################################################
# Performing the evaluation of the splsda feature selection algorithm
######################################################################################################



#gets the motif start location from the actual motifs table (TODO: Consider making it global if needed)
getMotifsStartLocations<-function(data_Motif_Locations)
{
  start_locations_list=c()
  for(i in 1:ncol(data_Motif_Locations))
  {
    start_locations_list= append(start_locations_list,data_Motif_Locations[3,i])
  }
  return(start_locations_list)
}

#gets the motif start location from the actual motifs table (TODO: Consider making it global if needed)
getMotifsEndLocations<-function(data_Motif_Locations)
{
  end_locations_list=c()
  for(i in 1:ncol(data_Motif_Locations))
  {
    end_locations_list= append(end_locations_list,data_Motif_Locations[4,i])
  }
  return(end_locations_list)
}


output_table_metrics=matrix(c("Run",rep(current_run_name,6)),nrow=7)
output_table_counts=matrix(c("Run",rep(current_run_name,5)),nrow=6)

output_table_metrics=cbind(output_table_metrics,c("Fold",rep(current_fold_name,6)))
output_table_counts=cbind(output_table_counts,c("Fold",rep(current_fold_name,5)))

#   output_table_metrics=matrix(c("Allele","all","all","all","all","all","all"),nrow=7)
#   output_table_counts=matrix(c("Allele","all","all","all","all","all"),nrow=6)


  output_table_metrics=cbind(output_table_metrics,c("Metric","sensitivity=TP/Real Positives=TP/(TP+FN)","specificity=TN/Real Negatives=TN/(TN+FP)",
            "precision=TP/splsda Positives=TP/(TP+FP)","neg_pred_value=TN/splsda Negatives=TN/(TN+FN)","J-index=Sensitivity + Specificity − 1","Accuracy"))
  output_table_counts=cbind(output_table_counts,c("Metric","TP","FP","TN","FN","Accuracy"))


#get the variable parameters column names
sample_dataset_names=output_file_names_short
print("sample_dataset_names")
print(sample_dataset_names)

col_names=sample_dataset_names


for(i in 1:length(actual_motifs_filenames))
{
  file_name=actual_motifs_filenames[i]
  data_Motif_Locations <- read.csv(file_name,header=TRUE,sep=",",colClasses="character")
  rownames(data_Motif_Locations)=NULL
  data_Motif_Locations=t(data_Motif_Locations)
  data_Motif_Locations=data_Motif_Locations[-1,]
  num_motifs=ncol(data_Motif_Locations)
  current_output_file_name_short=output_file_names_short[i]
  print("current_output_file_name_short")
  print(current_output_file_name_short)
  
  splsda_predicted_locations_file=paste(getSplsdaPredictedDiffLocations(training_dir),"/SPLSDA_Diff_locations_boolean_",current_output_file_name_short,".txt",sep="")
  splsda_predicted_locations=read.table(splsda_predicted_locations_file)
  splsda_predicted_locations=as.matrix(splsda_predicted_locations)
  positives_locations_all=which(splsda_predicted_locations==TRUE)  
  negatives_locations_all=which(splsda_predicted_locations==FALSE)
  total_positives_all=length(positives_locations_all)
  total_negatives_all=length(negatives_locations_all)
  
  splsda_testing_predicted_locations_file=paste(getSplsdaPredictedDiffLocations(testing_dir),"/SPLSDA_Diff_locations_boolean_",current_output_file_name_short,".txt",sep="")
  splsda_testing_predicted_locations=read.table(splsda_testing_predicted_locations_file)
  splsda_testing_predicted_locations=as.matrix(splsda_testing_predicted_locations)
  testing_positives_locations_all=which(splsda_testing_predicted_locations==TRUE)  
  testing_negatives_locations_all=which(splsda_testing_predicted_locations==FALSE)
  testing_total_positives_all=length(testing_positives_locations_all)
  testing_total_negatives_all=length(testing_negatives_locations_all)
  
if(OutputRuns)
{
  write(positives_locations_all,file=paste(getSplsdaPredictedDiffLocations(training_dir),"/SPLSDA_pos_locations_based_on_eval",current_output_file_name_short,".txt",sep=""),ncolumns = 1)
}

  
  
  TP_all=0
  FN_all=0
  
  motif_start_loc=getMotifsStartLocations(data_Motif_Locations)
  motif_end_loc=getMotifsEndLocations(data_Motif_Locations)
  
  for(j in 1:num_motifs)
  {
    curr_start_loc=motif_start_loc[j]
    curr_end_loc=motif_end_loc[j]
    curr_range=curr_start_loc:curr_end_loc
    
    current_motif_TP=length(which(positives_locations_all %in% curr_range))
    current_motif_TP_indexes=positives_locations_all[which(positives_locations_all%in%curr_range)]
    print("current_motif_TP_indexes")
    print(current_motif_TP_indexes)
    current_TP_loc_negative_in_testing=length(which(current_motif_TP_indexes%in%testing_negatives_locations_all))
    print("current_TP_loc_negative_in_testing")
    print(current_TP_loc_negative_in_testing)
    current_TP_loc_negative_in_testing_indexes=positives_locations_all[which(current_motif_TP_indexes%in%testing_negatives_locations_all)]
    print("current_TP_loc_negative_in_testing_indexes")
    print(current_TP_loc_negative_in_testing_indexes)    
    
    current_motif_TP=current_motif_TP-current_TP_loc_negative_in_testing
    print("current_motif_TP updated")
    print(current_motif_TP)
    total_positives_all=total_positives_all-current_TP_loc_negative_in_testing
    
    
    current_motif_FN=length(which(curr_range %in% negatives_locations_all))
    current_motif_FN=current_motif_FN+current_TP_loc_negative_in_testing
    print("current_motif_FN")
    print(current_motif_FN)
    total_negatives_all=total_negatives_all+current_TP_loc_negative_in_testing
    
    TP_all=TP_all+current_motif_TP
    FN_all=FN_all+current_motif_FN
  }
  

  FP_all=total_positives_all-TP_all
  TN_all=total_negatives_all-FN_all
  sensitivity_all=TP_all/(TP_all+FN_all)
  specificity_all=TN_all/(TN_all+FP_all)
  precision_all=TP_all/(TP_all+FP_all)
  negative_predictive_value_all=TN_all/(FN_all+TN_all)
  j_index_all=sensitivity_all+specificity_all-1  
  
  
  total=TP_all+TN_all+FP_all+FN_all
  accuracy_all=(TP_all+TN_all)/total
  
    output_table_metrics=cbind(output_table_metrics,c(col_names[i],round(sensitivity_all,digits=5),round(specificity_all,digits=5),
                                      round(precision_all,digits=5),round(negative_predictive_value_all,digits=5),round(j_index_all,digits=5),
                                      round(accuracy_all,digits=5)))
    output_table_counts=cbind(output_table_counts,c(col_names[i],TP_all,FP_all,TN_all,FN_all,accuracy_all))
}#End of for loop for every file names


#################################################
# Outputting final evaluation tables
#################################################

colnames(output_table_metrics)=output_table_metrics[1,]
output_table_metrics=output_table_metrics[-1,]
rownames(output_table_metrics)=NULL

colnames(output_table_counts)=output_table_counts[1,]
output_table_counts=output_table_counts[-1,]
rownames(output_table_counts)=NULL


if(OutputRuns)
{
  write.csv(output_table_metrics,file=splsda_differentiable_locations_rates)
  write.csv(output_table_counts,file=splsda_differentiable_locations_counts)
}


if(file.exists(evaluation_table_all_runs)==FALSE) {
  write.table(output_table_metrics,file=evaluation_table_all_runs,row.names=FALSE,append=TRUE,sep=",")
} else if(file.exists(evaluation_table_all_runs)==TRUE) {
  write.table(output_table_metrics,file=evaluation_table_all_runs,row.names=FALSE,append=TRUE,sep=",",col.names=FALSE)
}

if(file.exists(evaluation_counts_table_all_runs)==FALSE) {
  write.table(output_table_counts,file=evaluation_counts_table_all_runs,row.names=FALSE,append=TRUE,sep=",")
} else if(file.exists(evaluation_counts_table_all_runs)==TRUE) {
  write.table(output_table_counts,file=evaluation_counts_table_all_runs,row.names=FALSE,append=TRUE,sep=",",col.names=FALSE)
}


print("SplsdaEvaluation CODE ENDED!")   #print out in terminal an indication that the code has ended 

time=proc.time() - ptm  
print(paste("TIME in sec: ",time['elapsed']))

