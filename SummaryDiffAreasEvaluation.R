###############################################################
# Generating Average tables against all runs
###############################################################

print("Inside SummaryDiffAreasEvaluation.R")

OutputText=getOutputText()
OutputPlots=getOutputPlots()
param_list=getParamsForPlotTitle()
num_folds=getNumFolds()

print("Generating Differentiable Locations evaluations on all runs...")



for(dt in 1:length(data_type_vector))
{
  data_type=data_type_vector[dt]
  print("data_type")
  print(data_type)
  data_type_dir=""
  if(data_type=="Training")
  {
    data_type_dir=getTrainingAllOutputDir()
  }else{
    data_type_dir=getTestingAllOutputDir()
  }
  
  print("data_type_dir")
  print(data_type_dir)
  
  num_alleles=1 #ALL
  if(OutputPlots)
  {
    num_alleles=4 #0,1,2,ALL
  }
  num_metrics=6 #Spec, Sen, Acc, Prec, Recall,j-index
  num_stats=4 #Avg Max Min Stdiv
  num_runs=getNumRuns()
  #For differentiable locations
  all_evaluations_dir=paste(data_type_dir,"/Table_predicted_differentiable_locations_per_allele_metrics_all_runs.csv",sep="")
  all_evaluations_table=read.csv(all_evaluations_dir,header=FALSE,sep=",",stringsAsFactors=FALSE)
  
  print("all_evaluations_dir")
  print(all_evaluations_dir)

  
  
  all_evaluations_table=as.matrix(all_evaluations_table)
  colnames(all_evaluations_table)=all_evaluations_table[1,]
  all_evaluations_table=all_evaluations_table[-1,]
  
  print("all_evaluations_table")
  print(head(all_evaluations_table))
  
  average_table=matrix(nrow=num_alleles*num_metrics*num_stats)
  first_col=c(rep("Average",num_alleles*num_metrics),rep("Maximum",num_alleles*num_metrics),rep("Minimum",num_alleles*num_metrics),rep("Stand.Dev",num_alleles*num_metrics))
  average_table=cbind(average_table,first_col)
  second_2_columns=c(all_evaluations_table[1:(num_alleles*num_metrics),3],all_evaluations_table[1:(num_alleles*num_metrics),3],
                         all_evaluations_table[1:(num_alleles*num_metrics),3],all_evaluations_table[1:(num_alleles*num_metrics),3])
  print("length(second_2_columns)")
  print(length(second_2_columns))
  average_table=cbind(average_table,second_2_columns)
  average_table=average_table[,-1]
  row.names(average_table)=NULL
  
  if(OutputPlots)
  {
    allele0_SEN_indexes=seq(from=1,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele0_SPC_indexes=seq(from=2,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele0_PRC_indexes=seq(from=3,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele0_NPV_indexes=seq(from=4,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele0_JI_indexes=seq(from=5,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele0_ACC_indexes=seq(from=6,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    
    allele1_SEN_indexes=seq(from=7,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele1_SPC_indexes=seq(from=8,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele1_PRC_indexes=seq(from=9,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele1_NPV_indexes=seq(from=10,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele1_JI_indexes=seq(from=11,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele1_ACC_indexes=seq(from=12,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    
    allele2_SEN_indexes=seq(from=13,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele2_SPC_indexes=seq(from=14,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele2_PRC_indexes=seq(from=15,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele2_NPV_indexes=seq(from=16,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele2_JI_indexes=seq(from=17,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    allele2_ACC_indexes=seq(from=18,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    
    all_SEN_indexes=seq(from=19,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_SPC_indexes=seq(from=20,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_PRC_indexes=seq(from=21,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_NPV_indexes=seq(from=22,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_JI_indexes=seq(from=23,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_ACC_indexes=seq(from=24,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    
    
    for(c in 3:ncol(all_evaluations_table))
    {
      curr_column=all_evaluations_table[,c]
      #   Computing Averages
      avg_allele0_SEN=mean(as.numeric(curr_column[allele0_SEN_indexes]),na.rm=TRUE)
      avg_allele0_SPC=mean(as.numeric(curr_column[allele0_SPC_indexes]),na.rm=TRUE)
      avg_allele0_PRC=mean(as.numeric(curr_column[allele0_PRC_indexes]),na.rm=TRUE)
      avg_allele0_NPV=mean(as.numeric(curr_column[allele0_NPV_indexes]),na.rm=TRUE)
      avg_allele0_JI=mean(as.numeric(curr_column[allele0_JI_indexes]),na.rm=TRUE)
      avg_allele0_ACC=mean(as.numeric(curr_column[allele0_ACC_indexes]),na.rm=TRUE)
      
      avg_allele1_SEN=mean(as.numeric(curr_column[allele1_SEN_indexes]),na.rm=TRUE)
      avg_allele1_SPC=mean(as.numeric(curr_column[allele1_SPC_indexes]),na.rm=TRUE)
      avg_allele1_PRC=mean(as.numeric(curr_column[allele1_PRC_indexes]),na.rm=TRUE)
      avg_allele1_NPV=mean(as.numeric(curr_column[allele1_NPV_indexes]),na.rm=TRUE)
      avg_allele1_JI=mean(as.numeric(curr_column[allele1_JI_indexes]),na.rm=TRUE)
      avg_allele1_ACC=mean(as.numeric(curr_column[allele1_ACC_indexes]),na.rm=TRUE)
      
      avg_allele2_SEN=mean(as.numeric(curr_column[allele2_SEN_indexes]),na.rm=TRUE)
      avg_allele2_SPC=mean(as.numeric(curr_column[allele2_SPC_indexes]),na.rm=TRUE)
      avg_allele2_PRC=mean(as.numeric(curr_column[allele2_PRC_indexes]),na.rm=TRUE)
      avg_allele2_NPV=mean(as.numeric(curr_column[allele2_NPV_indexes]),na.rm=TRUE)
      avg_allele2_JI=mean(as.numeric(curr_column[allele2_JI_indexes]),na.rm=TRUE)
      avg_allele2_ACC=mean(as.numeric(curr_column[allele2_ACC_indexes]),na.rm=TRUE)
      
      avg_all_SEN=mean(as.numeric(curr_column[all_SEN_indexes]),na.rm=TRUE)
      avg_all_SPC=mean(as.numeric(curr_column[all_SPC_indexes]),na.rm=TRUE)
      avg_all_PRC=mean(as.numeric(curr_column[all_PRC_indexes]),na.rm=TRUE)
      avg_all_NPV=mean(as.numeric(curr_column[all_NPV_indexes]),na.rm=TRUE)
      avg_all_JI=mean(as.numeric(curr_column[all_JI_indexes]),na.rm=TRUE)
      avg_all_ACC=mean(as.numeric(curr_column[all_ACC_indexes]),na.rm=TRUE)
      
      
      
      #   Computing Maximum
      max_allele0_SEN=as.numeric(max(curr_column[allele0_SEN_indexes],,na.rm=TRUE))
      max_allele0_SPC=as.numeric(max(curr_column[allele0_SPC_indexes],na.rm=TRUE))
      max_allele0_PRC=as.numeric(max(curr_column[allele0_PRC_indexes],na.rm=TRUE))
      max_allele0_NPV=as.numeric(max(curr_column[allele0_NPV_indexes],na.rm=TRUE))
      max_allele0_JI=as.numeric(max(curr_column[allele0_JI_indexes],na.rm=TRUE))
      max_allele0_ACC=as.numeric(max(curr_column[allele0_ACC_indexes],na.rm=TRUE))
      
      max_allele1_SEN=as.numeric(max(curr_column[allele1_SEN_indexes],na.rm=TRUE))
      max_allele1_SPC=as.numeric(max(curr_column[allele1_SPC_indexes],na.rm=TRUE))
      max_allele1_PRC=as.numeric(max(curr_column[allele1_PRC_indexes],na.rm=TRUE))
      max_allele1_NPV=as.numeric(max(curr_column[allele1_NPV_indexes],na.rm=TRUE))
      max_allele1_JI=as.numeric(max(curr_column[allele1_JI_indexes],na.rm=TRUE))
      max_allele1_ACC=as.numeric(max(curr_column[allele1_ACC_indexes],na.rm=TRUE))
      
      max_allele2_SEN=as.numeric(max(curr_column[allele2_SEN_indexes],na.rm=TRUE))
      max_allele2_SPC=as.numeric(max(curr_column[allele2_SPC_indexes],na.rm=TRUE))
      max_allele2_PRC=as.numeric(max(curr_column[allele2_PRC_indexes],na.rm=TRUE))
      max_allele2_NPV=as.numeric(max(curr_column[allele2_NPV_indexes],na.rm=TRUE))
      max_allele2_JI=as.numeric(max(curr_column[allele2_JI_indexes],na.rm=TRUE))
      max_allele2_ACC=as.numeric(max(curr_column[allele2_ACC_indexes],na.rm=TRUE))
      
      max_all_SEN=as.numeric(max(curr_column[all_SEN_indexes],na.rm=TRUE))
      max_all_SPC=as.numeric(max(curr_column[all_SPC_indexes],na.rm=TRUE))
      max_all_PRC=as.numeric(max(curr_column[all_PRC_indexes],na.rm=TRUE))
      max_all_NPV=as.numeric(max(curr_column[all_NPV_indexes],na.rm=TRUE))
      max_all_JI=as.numeric(max(curr_column[all_JI_indexes],na.rm=TRUE))
      max_all_ACC=as.numeric(max(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      #   Computing Minimum
      min_allele0_SEN=as.numeric(min(curr_column[allele0_SEN_indexes],na.rm=TRUE))
      min_allele0_SPC=as.numeric(min(curr_column[allele0_SPC_indexes],na.rm=TRUE))
      min_allele0_PRC=as.numeric(min(curr_column[allele0_PRC_indexes],na.rm=TRUE))
      min_allele0_NPV=as.numeric(min(curr_column[allele0_NPV_indexes],na.rm=TRUE))
      min_allele0_JI=as.numeric(min(curr_column[allele0_JI_indexes],na.rm=TRUE))
      min_allele0_ACC=as.numeric(min(curr_column[allele0_ACC_indexes],na.rm=TRUE))
      
      min_allele1_SEN=as.numeric(min(curr_column[allele1_SEN_indexes],na.rm=TRUE))
      min_allele1_SPC=as.numeric(min(curr_column[allele1_SPC_indexes],na.rm=TRUE))
      min_allele1_PRC=as.numeric(min(curr_column[allele1_PRC_indexes],na.rm=TRUE))
      min_allele1_NPV=as.numeric(min(curr_column[allele1_NPV_indexes],na.rm=TRUE))
      min_allele1_JI=as.numeric(min(curr_column[allele1_JI_indexes],na.rm=TRUE))
      min_allele1_ACC=as.numeric(min(curr_column[allele1_ACC_indexes],na.rm=TRUE))
      
      min_allele2_SEN=as.numeric(min(curr_column[allele2_SEN_indexes],na.rm=TRUE))
      min_allele2_SPC=as.numeric(min(curr_column[allele2_SPC_indexes],na.rm=TRUE))
      min_allele2_PRC=as.numeric(min(curr_column[allele2_PRC_indexes],na.rm=TRUE))
      min_allele2_NPV=as.numeric(min(curr_column[allele2_NPV_indexes],na.rm=TRUE))
      min_allele2_JI=as.numeric(min(curr_column[allele2_JI_indexes],na.rm=TRUE))
      min_allele2_ACC=as.numeric(min(curr_column[allele2_ACC_indexes],na.rm=TRUE))
      
      min_all_SEN=as.numeric(min(curr_column[all_SEN_indexes],na.rm=TRUE))
      min_all_SPC=as.numeric(min(curr_column[all_SPC_indexes],na.rm=TRUE))
      min_all_PRC=as.numeric(min(curr_column[all_PRC_indexes],na.rm=TRUE))
      min_all_NPV=as.numeric(min(curr_column[all_NPV_indexes],na.rm=TRUE))
      min_all_JI=as.numeric(min(curr_column[all_JI_indexes],na.rm=TRUE))
      min_all_ACC=as.numeric(min(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      #   Computing Standard deviation
      
      sd_allele0_SEN=as.numeric(sd(curr_column[allele0_SEN_indexes],na.rm=TRUE))
      sd_allele0_SPC=as.numeric(sd(curr_column[allele0_SPC_indexes],na.rm=TRUE))
      sd_allele0_PRC=as.numeric(sd(curr_column[allele0_PRC_indexes],na.rm=TRUE))
      sd_allele0_NPV=as.numeric(sd(curr_column[allele0_NPV_indexes],na.rm=TRUE))
      sd_allele0_JI=as.numeric(sd(curr_column[allele0_JI_indexes],na.rm=TRUE))
      sd_allele0_ACC=as.numeric(sd(curr_column[allele0_ACC_indexes],na.rm=TRUE))
      
      sd_allele1_SEN=as.numeric(sd(curr_column[allele1_SEN_indexes],na.rm=TRUE))
      sd_allele1_SPC=as.numeric(sd(curr_column[allele1_SPC_indexes],na.rm=TRUE))
      sd_allele1_PRC=as.numeric(sd(curr_column[allele1_PRC_indexes],na.rm=TRUE))
      sd_allele1_NPV=as.numeric(sd(curr_column[allele1_NPV_indexes],na.rm=TRUE))
      sd_allele1_JI=as.numeric(sd(curr_column[allele1_JI_indexes],na.rm=TRUE))
      sd_allele1_ACC=as.numeric(sd(curr_column[allele1_ACC_indexes],na.rm=TRUE))
      
      sd_allele2_SEN=as.numeric(sd(curr_column[allele2_SEN_indexes],na.rm=TRUE))
      sd_allele2_SPC=as.numeric(sd(curr_column[allele2_SPC_indexes],na.rm=TRUE))
      sd_allele2_PRC=as.numeric(sd(curr_column[allele2_PRC_indexes],na.rm=TRUE))
      sd_allele2_NPV=as.numeric(sd(curr_column[allele2_NPV_indexes],na.rm=TRUE))
      sd_allele2_JI=as.numeric(sd(curr_column[allele2_JI_indexes],na.rm=TRUE))
      sd_allele2_ACC=as.numeric(sd(curr_column[allele2_ACC_indexes],na.rm=TRUE))
      
      sd_all_SEN=as.numeric(sd(curr_column[all_SEN_indexes],na.rm=TRUE))
      sd_all_SPC=as.numeric(sd(curr_column[all_SPC_indexes],na.rm=TRUE))
      sd_all_PRC=as.numeric(sd(curr_column[all_PRC_indexes],na.rm=TRUE))
      sd_all_NPV=as.numeric(sd(curr_column[all_NPV_indexes],na.rm=TRUE))
      sd_all_JI=as.numeric(sd(curr_column[all_JI_indexes],na.rm=TRUE))
      sd_all_ACC=as.numeric(sd(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      
      average_curr_column=c(round(avg_allele0_SEN,digits=5),round(avg_allele0_SPC,digits=5),round(avg_allele0_PRC,digits=5)
                            ,round(avg_allele0_NPV,digits=5),round(avg_allele0_JI,digits=5),round(avg_allele0_ACC,digits=5),round(avg_allele1_SEN,digits=5)
                            ,round(avg_allele1_SPC,digits=5),round(avg_allele1_PRC,digits=5),round(avg_allele1_NPV,digits=5)
                            ,round(avg_allele1_JI,digits=5),round(avg_allele1_ACC,digits=5),round(avg_allele2_SEN,digits=5),round(avg_allele2_SPC,digits=5)
                            ,round(avg_allele2_PRC,digits=5),round(avg_allele2_NPV,digits=5),round(avg_allele2_JI,digits=5),round(avg_allele2_ACC,digits=5)
                            ,round(avg_all_SEN,digits=5),round(avg_all_SPC,digits=5),round(avg_all_PRC,digits=5)
                            ,round(avg_all_NPV,digits=5),round(avg_all_JI,digits=5),round(avg_all_ACC,digits=5)
                            ,round(max_allele0_SEN,digits=5),round(max_allele0_SPC,digits=5),round(max_allele0_PRC,digits=5)
                            ,round(max_allele0_NPV,digits=5),round(max_allele0_JI,digits=5),round(max_allele0_ACC,digits=5),round(max_allele1_SEN,digits=5)
                            ,round(max_allele1_SPC,digits=5),round(max_allele1_PRC,digits=5),round(max_allele1_NPV,digits=5),round(max_allele1_JI,digits=5)
                            ,round(max_allele1_ACC,digits=5),round(max_allele2_SEN,digits=5),round(max_allele2_SPC,digits=5)
                            ,round(max_allele2_PRC,digits=5),round(max_allele2_NPV,digits=5),round(max_allele2_JI,digits=5),round(max_allele2_ACC,digits=5)
                            ,round(max_all_SEN,digits=5),round(max_all_SPC,digits=5)
                            ,round(max_all_PRC,digits=5),round(max_all_NPV,digits=5),round(max_all_JI,digits=5),round(max_all_ACC,digits=5)
                            ,round(min_allele0_SEN,digits=5),round(min_allele0_SPC,digits=5),round(min_allele0_PRC,digits=5)
                            ,round(min_allele0_NPV,digits=5),round(min_allele0_JI,digits=5),round(min_allele0_ACC,digits=5),round(min_allele1_SEN,digits=5)
                            ,round(min_allele1_SPC,digits=5),round(min_allele1_PRC,digits=5),round(min_allele1_NPV,digits=5),round(min_allele1_JI,digits=5)
                            ,round(min_allele1_ACC,digits=5),round(min_allele2_SEN,digits=5),round(min_allele2_SPC,digits=5)
                            ,round(min_allele2_PRC,digits=5),round(min_allele2_NPV,digits=5),round(min_allele2_ACC,digits=5),round(min_allele2_JI,digits=5)
                            ,round(min_all_SEN,digits=5),round(min_all_SPC,digits=5)
                            ,round(min_all_PRC,digits=5),round(min_all_NPV,digits=5),round(min_all_JI,digits=5),round(min_all_ACC,digits=5)
                            ,round(sd_allele0_SEN,digits=5),round(sd_allele0_SPC,digits=5),round(sd_allele0_PRC,digits=5)
                            ,round(sd_allele0_NPV,digits=5),round(sd_allele0_JI,digits=5),round(sd_allele0_ACC,digits=5),round(sd_allele1_SEN,digits=5)
                            ,round(sd_allele1_SPC,digits=5),round(sd_allele1_PRC,digits=5),round(sd_allele1_NPV,digits=5),round(sd_allele1_JI,digits=5)
                            ,round(sd_allele1_ACC,digits=5),round(sd_allele2_SEN,digits=5),round(sd_allele2_SPC,digits=5)
                            ,round(sd_allele2_PRC,digits=5),round(sd_allele2_NPV,digits=5),round(sd_allele2_JI,digits=5),round(sd_allele2_ACC,digits=5)
                            ,round(sd_all_SEN,digits=5),round(sd_all_SPC,digits=5)
                            ,round(sd_all_PRC,digits=5),round(sd_all_NPV,digits=5),round(sd_all_JI,digits=5),round(sd_all_ACC,digits=5))
      
      average_table=cbind(average_table,average_curr_column)
      
    }
  }else{
    
    all_SEN_indexes=seq(from=1,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_SPC_indexes=seq(from=2,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_PRC_indexes=seq(from=3,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_NPV_indexes=seq(from=4,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_JI_indexes=seq(from=5,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    all_ACC_indexes=seq(from=6,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
    
    
    for(c in 4:ncol(all_evaluations_table))
    {
      curr_column=all_evaluations_table[,c]
      #   Computing Averages
      
      avg_all_SEN=mean(as.numeric(curr_column[all_SEN_indexes]),na.rm=TRUE)
      avg_all_SPC=mean(as.numeric(curr_column[all_SPC_indexes]),na.rm=TRUE)
      avg_all_PRC=mean(as.numeric(curr_column[all_PRC_indexes]),na.rm=TRUE)
      avg_all_NPV=mean(as.numeric(curr_column[all_NPV_indexes]),na.rm=TRUE)
      avg_all_JI=mean(as.numeric(curr_column[all_JI_indexes]),na.rm=TRUE)
      avg_all_ACC=mean(as.numeric(curr_column[all_ACC_indexes]),na.rm=TRUE)
      
      #   Computing Maximum
      
      max_all_SEN=as.numeric(max(curr_column[all_SEN_indexes],na.rm=TRUE))
      max_all_SPC=as.numeric(max(curr_column[all_SPC_indexes],na.rm=TRUE))
      max_all_PRC=as.numeric(max(curr_column[all_PRC_indexes],na.rm=TRUE))
      max_all_NPV=as.numeric(max(curr_column[all_NPV_indexes],na.rm=TRUE))
      max_all_JI=as.numeric(max(curr_column[all_JI_indexes],na.rm=TRUE))
      max_all_ACC=as.numeric(max(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      #   Computing Minimum
      
      min_all_SEN=as.numeric(min(curr_column[all_SEN_indexes],na.rm=TRUE))
      min_all_SPC=as.numeric(min(curr_column[all_SPC_indexes],na.rm=TRUE))
      min_all_PRC=as.numeric(min(curr_column[all_PRC_indexes],na.rm=TRUE))
      min_all_NPV=as.numeric(min(curr_column[all_NPV_indexes],na.rm=TRUE))
      min_all_JI=as.numeric(min(curr_column[all_JI_indexes],na.rm=TRUE))
      min_all_ACC=as.numeric(min(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      #   Computing Standard deviation
      
      sd_all_SEN=as.numeric(sd(curr_column[all_SEN_indexes],na.rm=TRUE))
      sd_all_SPC=as.numeric(sd(curr_column[all_SPC_indexes],na.rm=TRUE))
      sd_all_PRC=as.numeric(sd(curr_column[all_PRC_indexes],na.rm=TRUE))
      sd_all_NPV=as.numeric(sd(curr_column[all_NPV_indexes],na.rm=TRUE))
      sd_all_JI=as.numeric(sd(curr_column[all_JI_indexes],na.rm=TRUE))
      sd_all_ACC=as.numeric(sd(curr_column[all_ACC_indexes],na.rm=TRUE))
      
      
      
      average_curr_column=c(round(avg_all_SEN,digits=5),round(avg_all_SPC,digits=5),round(avg_all_PRC,digits=5)
                            ,round(avg_all_NPV,digits=5),round(avg_all_JI,digits=5),round(avg_all_ACC,digits=5)
                            ,round(max_all_SEN,digits=5),round(max_all_SPC,digits=5)
                            ,round(max_all_PRC,digits=5),round(max_all_NPV,digits=5),round(max_all_JI,digits=5)
                            ,round(max_all_ACC,digits=5),round(min_all_SEN,digits=5),round(min_all_SPC,digits=5)
                            ,round(min_all_PRC,digits=5),round(min_all_NPV,digits=5),round(min_all_JI,digits=5)
                            ,round(min_all_ACC,digits=5),round(sd_all_SEN,digits=5),round(sd_all_SPC,digits=5)
                            ,round(sd_all_PRC,digits=5),round(sd_all_NPV,digits=5),round(sd_all_JI,digits=5)
                            ,round(sd_all_ACC,digits=5))
      
      average_table=cbind(average_table,average_curr_column)
    }
  }
  all_evaluations_table_colnames=colnames(all_evaluations_table)
  colnames(average_table)=c(" ",all_evaluations_table_colnames[3:length(all_evaluations_table_colnames)])
  
  print("avg_output_file")
  print(avg_output_file)
  print("head(average_table)")
  print(head(average_table))
  
  avg_output_file=paste(data_type_dir,"/Table_predicted_differentiable_locations_per_allele_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
  write.csv(average_table,file=avg_output_file)
  
}#end of for(dt in 1:length(data_type_vector))






#################################################################
# Reformatting Tables to plot graphs
#################################################################

print("Changing format of table and plotting evaluation of differentiable locations prediction on all runs")

#GET TRAINING TABLES

data_type="Training"
data_type_dir=getTrainingAllOutputDir()

splsda_avg_output_file_train=paste(data_type_dir,"/Table_splsda_predicted_differentiable_locations_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
splsda_avg_table_train=read.csv(splsda_avg_output_file_train,header=FALSE,sep=",",stringsAsFactors=FALSE)

print("splsda_avg_output_file_train")
print(splsda_avg_output_file_train)

splsda_avg_table_train=as.matrix(splsda_avg_table_train)
colnames(splsda_avg_table_train)=splsda_avg_table_train[1,]
splsda_avg_table_train=splsda_avg_table_train[-1,]
splsda_avg_table_train=splsda_avg_table_train[,-1]

print("head(splsda_avg_table_train)")
print(head(splsda_avg_table_train))



average_table_train_file=paste(data_type_dir,"/Table_predicted_differentiable_locations_per_allele_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
average_table_train=read.csv(average_table_train_file,header=FALSE,sep=",",stringsAsFactors=FALSE)

print("average_table_train_file")
print(average_table_train_file)

average_table_train=as.matrix(average_table_train)
colnames(average_table_train)=average_table_train[1,]
average_table_train=average_table_train[-1,]
average_table_train=average_table_train[,-1]

print("head(average_table_train)")
print(head(average_table_train))






#GET TESTING TABLES

data_type="testing"
data_type_dir=getTestingAllOutputDir()

splsda_avg_output_file_test=paste(data_type_dir,"/Table_splsda_predicted_differentiable_locations_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
splsda_avg_table_test=read.csv(splsda_avg_output_file_test,header=FALSE,sep=",",stringsAsFactors=FALSE)

print("splsda_avg_output_file_test")
print(splsda_avg_output_file_test)

splsda_avg_table_test=as.matrix(splsda_avg_table_test)
colnames(splsda_avg_table_test)=splsda_avg_table_test[1,]
splsda_avg_table_test=splsda_avg_table_test[-1,]
splsda_avg_table_test=splsda_avg_table_test[,-1]

print("head(splsda_avg_table_test)")
print(head(splsda_avg_table_test))


average_table_test_file=paste(data_type_dir,"/Table_predicted_differentiable_locations_per_allele_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
average_table_test=read.csv(average_table_test_file,header=FALSE,sep=",",stringsAsFactors=FALSE)

print("average_table_test_file")
print(average_table_test_file)

average_table_test=as.matrix(average_table_test)
colnames(average_table_test)=average_table_test[1,]
average_table_test=average_table_test[-1,]
average_table_test=average_table_test[,-1]

print("head(average_table_test)")
print(head(average_table_test))














variables=colnames(average_table_test)[-1:-2]
variables=variables[-length(variables)]
variables=unlist(strsplit(variables,"_",fixed=TRUE))
variable_name=variables[1]
index_to_keep=seq(from=0,to=length(variables),by=2)
variables=variables[index_to_keep]
num_variables=length(variables)
if(OutputText)
{
  variables=rep(variables,times=4)
}

graph_table=data.frame()

num_metrics=6

all_stdv_col=as.numeric(c(average_table_train[(1+18),3:(num_variables+2)],splsda_avg_table_train[(1+18),3:(num_variables+2)],
                    average_table_train[(2+18),3:(num_variables+2)],splsda_avg_table_train[(2+18),3:(num_variables+2)],
                    average_table_train[(3+18),3:(num_variables+2)],splsda_avg_table_train[(3+18),3:(num_variables+2)],
                    average_table_train[(4+18),3:(num_variables+2)],splsda_avg_table_train[(4+18),3:(num_variables+2)],
                    average_table_train[(5+18),3:(num_variables+2)],splsda_avg_table_train[(5+18),3:(num_variables+2)],
                    average_table_train[(6+18),3:(num_variables+2)],splsda_avg_table_train[(6+18),3:(num_variables+2)],
                    average_table_test[(1+18),3:(num_variables+2)],splsda_avg_table_test[(1+18),3:(num_variables+2)],
                    average_table_test[(2+18),3:(num_variables+2)],splsda_avg_table_test[(2+18),3:(num_variables+2)],
                    average_table_test[(3+18),3:(num_variables+2)],splsda_avg_table_test[(3+18),3:(num_variables+2)],
                    average_table_test[(4+18),3:(num_variables+2)],splsda_avg_table_test[(4+18),3:(num_variables+2)],
                    average_table_test[(5+18),3:(num_variables+2)],splsda_avg_table_test[(5+18),3:(num_variables+2)],
                    average_table_test[(6+18),3:(num_variables+2)],splsda_avg_table_test[(6+18),3:(num_variables+2)]))
print("all_stdv_col")
print(all_stdv_col)

max_stdv=as.numeric(max(all_stdv_col,na.rm=TRUE))
print("max_stdv")
print(max_stdv)
y_lim=1+max_stdv
print("y_lim")
print(y_lim)

if(OutputPlots)
{
  for( i in 1:num_metrics)
  {
    genotype_col=c(rep(0,times=num_variables),rep(1,times=num_variables),rep(2,times=num_variables),rep("all",times=num_variables))    
    variable_col=variables
    variable_col=as.numeric(variable_col)
    metric=unlist(strsplit(as.character(average_table[i,2]),"="))[1]
    avg_col=as.numeric(c(average_table[i,3:(num_variables+2)],average_table[(i+5),3:(num_variables+2)],average_table[(i+10),3:(num_variables+2)],average_table[(i+15),3:(num_variables+2)]))
    max_col=as.numeric(c(average_table[(i+20),3:(num_variables+2)],average_table[(i+25),3:(num_variables+2)],average_table[(i+30),3:(num_variables+2)],average_table[(i+35),3:(num_variables+2)]))
    min_col=as.numeric(c(average_table[(i+40),3:(num_variables+2)],average_table[(i+45),3:(num_variables+2)],average_table[(i+50),3:(num_variables+2)],average_table[(i+55),3:(num_variables+2)]))
    stdv_col=as.numeric(c(average_table[(i+60),3:(num_variables+2)],average_table[(i+65),3:(num_variables+2)],average_table[(i+70),3:(num_variables+2)],average_table[(i+75),3:(num_variables+2)]))
    
    graph_table=data.frame(Genotype=factor(genotype_col),Variable_Parameter=factor(variable_col),Average=avg_col,Maximum=max_col,Minimum=min_col,
                           Standard_Deviation=stdv_col)
    col_names=c("Genotype","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
    colnames(graph_table)=col_names
    Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
    Genotype=graph_table$Genotype
    error_bar_min=graph_table$Average-graph_table$Standard_Deviation
    error_bar_min[which(error_bar_min<0)]=0
    error_bar_max=graph_table$Average+graph_table$Standard_Deviation
    error_bar_max[which(error_bar_max>1)]=1
    plots_dir=paste(data_type_dir,"/Graph_predicted_differentiable_locations_per_allele_",metric,".jpeg",sep="")
    plot=qplot(graph_table$Variable_Parameter,graph_table$Average, data=graph_table, colour=Genotype)
    plot=plot+geom_line(aes(group=graph_table$Genotype))
    plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
    plot=plot+ggtitle(paste(metric,"of Differentiable SNP Locations\nVS.",variable_name,"Over",num_runs,"Runs\n",param_list,data_type))
    plot=plot+theme(plot.title = element_text(size = 15))  
    plot=plot+labs(x=variable_name,y=paste("Average",metric))
    plot=plot+ scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
    plot=plot+scale_y_continuous(limit=c(0,1))
    ggsave(plots_dir)
  }
}else {
  for( i in 1:num_metrics)
  {
    data_type_col=c(rep("Training",times=num_variables*2),rep("Testing",times=num_variables*2))
    print("data_type_col")
    print(data_type_col)
    alg_col=c(rep("SNP-DMF",times=num_variables),rep("SPLSDA",times=num_variables))
    alg_col=rep(alg_col,2)# once for training and once for testing
    print("alg_col")
    print(alg_col)
    
    data_type_alg_col=apply(cbind(data_type_col, alg_col), 1, function(x) paste(x, collapse="_"))
    print("data_type_alg_col")
    print(data_type_alg_col)
    
    variable_col=rep(variables,4)
    variable_col=as.numeric(variable_col)
    print("variable_col")
    print(variable_col)
    metric=unlist(strsplit(as.character(average_table_test[i,2]),"="))[1]
    avg_col=as.numeric(c(average_table_train[i,3:(num_variables+2)],splsda_avg_table_train[i,3:(num_variables+2)],
                         average_table_test[i,3:(num_variables+2)],splsda_avg_table_test[i,3:(num_variables+2)]))
    max_col=as.numeric(c(average_table_train[(i+6),3:(num_variables+2)],splsda_avg_table_train[(i+6),3:(num_variables+2)],
                         average_table_test[(i+6),3:(num_variables+2)],splsda_avg_table_test[(i+6),3:(num_variables+2)]))
    min_col=as.numeric(c(average_table_train[(i+12),3:(num_variables+2)],splsda_avg_table_train[(i+12),3:(num_variables+2)],
                         average_table_test[(i+12),3:(num_variables+2)],splsda_avg_table_test[(i+12),3:(num_variables+2)]))
    stdv_col=as.numeric(c(average_table_train[(i+18),3:(num_variables+2)],splsda_avg_table_train[(i+18),3:(num_variables+2)],
                          average_table_test[(i+18),3:(num_variables+2)],splsda_avg_table_test[(i+18),3:(num_variables+2)]))

    print("avg_col")
    print(avg_col)
    print("max_col")
    print(max_col)
    print("min_col")
    print(min_col)
    print("stdv_col")
    print(stdv_col)
    
    print("y_lim")
    print(y_lim)
#     graph_table=data.frame(Data_Type=factor(data_type_col),Algorithm=factor(alg_col),Variable_Parameter=factor(variable_col),Average=avg_col,Maximum=max_col,Minimum=min_col,
#                            Standard_Deviation=stdv_col)
    graph_table=data.frame(Data_Type_Alg=factor(data_type_alg_col),Data_Type=factor(data_type_col),Algorithm=factor(alg_col),Variable_Parameter=factor(variable_col),Average=avg_col,Maximum=max_col,Minimum=min_col,
                           Standard_Deviation=stdv_col)
    col_names=c("Data_Type_Alg","Data_Type","Algorithm","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
    colnames(graph_table)=col_names
#     graph_table_file=paste(getTrainingTestingPerformanceEvaluationDir(),"/predicted_differentiable_locations_benchmarks_metrics_all_runs_MAX_MIN_AVG_STDV_Reformatted",metric,".csv",sep="")
#     write.csv(graph_table,file=graph_table_file)

#################################################################
# Plotting Graphs for diff locations
#################################################################

    Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
print("Variable_Parameter")
print(Variable_Parameter)
    Data_Type_Alg=graph_table$Data_Type_Alg
print("Data_Type_Alg")
print(Data_Type_Alg)
    Data_Type=graph_table$Data_Type
print("Data_Type")
print(Data_Type)
    Algorithm=graph_table$Algorithm
print("Algorithm")
print(Algorithm)
    error_bar_min=graph_table$Average-graph_table$Standard_Deviation
    error_bar_min[which(error_bar_min<0)]=0
    error_bar_max=graph_table$Average+graph_table$Standard_Deviation
    plots_dir=paste(getTrainingTestingPerformanceEvaluationDir(),"/Graph_predicted_differentiable_locations_benchmarks_",metric,".jpeg",sep="")
    plot=qplot(graph_table$Variable_Parameter,graph_table$Average, data=graph_table, colour=Algorithm)
    plot=plot+geom_line(aes(group=graph_table$Data_Type_Alg,linetype=Data_Type))    
# plot=qplot(graph_table$Variable_Parameter,graph_table$Average, data=graph_table, colour=Algorithm)
#     plot=plot+geom_line(aes(group=graph_table$Data_Type_Alg,linetype=graph_table$Data_Type))
    plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
    plot=plot+ggtitle(paste(metric,"of Differentiable SNP Locations\nVS.",variable_name,"Over",num_runs,"Runs","and ",num_folds,"Folds\n",param_list))
    plot=plot+theme(plot.title = element_text(size = 15))  
    plot=plot+labs(x=variable_name,y=paste("Average",metric))
    plot=plot+ scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
    plot=plot+scale_y_continuous(limit=c(0,y_lim))
    ggsave(plots_dir)
  }
}






