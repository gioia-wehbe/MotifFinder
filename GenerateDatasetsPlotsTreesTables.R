#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# GenerateDatasetsPlotsTreesTables: 
# A program that runs DataGenerator, followed by PlotGenerator on many different input files. 
# To run it on the different input files, set the working directory to the directory where GenerateDatasetsAndPlots.R is
# and run GenerateDatasetsAndPlots.R while passing the directory where all the runs files are located along with the variable
# parameter currently being worked on
# Author: Gioia Wehbe
# LAST UPDATED: 23-02-15
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print("Inside GenerateDatasetsPlotsTreesTables.R")

# ptm <- proc.time()
# args <- commandArgs(trailingOnly = TRUE)
# 
# datasets_dir=args[1]
# runs_dir=args[2]
# variable_parameter=args[3]      #name of the parameter that will be changing
# pmin_list=args[4]               #Probability of none-occurence parameter
# eval_shift_window=as.numeric(args[5])        #=1: number of SNPs allowed to be added to the ends of a predicted motifs so it is stil considered T (used in EvaluationFullDiffMotifs)
# eval_diss_threshold=as.numeric(args[6])     #=2: number of max allowed insertions, deletions and replacements when evaluation motifs (used in EvaluationFullDiffMotifs)
# fold_diff_upper_bound=as.numeric(args[7])   #=2
# fold_diff_lower_bound=as.numeric(args[8])   #=0.5
# min_motif_length=as.numeric(args[9])        #4
# 
# OutputText=args[10]              #if FALSE => do not output least needed items, if true => output everything
# OutputPlots=args[11]             #if FALSE => do not run PlotGenerator & Generate plots, if TRUE => plot them
# OutputRuns<<-args[12]            #if FALSE => do not output Runs Output, if true => output everything
# num_runs=as.numeric(args[13])
# num_folds=as.numeric(args[14])






#for TESTING!!!!
datasets_dir="C:\\Users\\gioia.wehbe\\Experiments Simulated Data\\Experiment\\ProbaNoneOcc\\SimulatedDataSets"
runs_dir="C:\\Users\\gioia.wehbe\\Experiments Simulated Data\\Experiment\\ProbaNoneOcc\\RunsOutput"
variable_parameter="ProbaNoneOcc"      #name of the parameter that will be changing
pmin_list="0.8"
eval_shift_window=1       #=1: number of SNPs allowed to be added to the ends of a predicted motifs so it is stil considered T (used in EvaluationFullDiffMotifs)
eval_diss_threshold=2     #=2: number of max allowed insertions, deletions and replacements when evaluation motifs (used in EvaluationFullDiffMotifs)
fold_diff_upper_bound=2
fold_diff_lower_bound=0.5
min_motif_length=4
OutputText="FALSE"     #if FALSE => do not output least needed items, if true => output everything
OutputPlots="FALSE"    #if FALSE => do not run PlotGenerator & Generate plots, if TRUE => plot them
OutputRuns<<-"FALSE"   #if FALSE => do not output Runs Output, if true => output everything
num_runs=10
num_folds=10
setwd("C:\\Users\\gioia.wehbe\\Rprojects\\ThesisSrcCombined(18-06-15)")




library("ggplot2")
library("TraMineR") 
library("PST") 
library("plyr")
library("caret")
library('ROCR')
library("mixOmics")

home_dir=getwd()

getNumRuns<<-function()
{
  return(num_runs)
}

getNumFolds<<-function()
{
  return(num_folds)
}

getHomeDir<<-function()
{
  return(home_dir)
}

getPminListText<<-function()
{
  return(pmin_list)
}


getPminList<<-function()
{
  pmin_list_num=unlist(strsplit(pmin_list,","))
  return(pmin_list_num)
}

getEvaluationShiftWindow<<-function()
{
  return(eval_shift_window)
}

getEvaluationDistanceThreshold<<-function()
{
  return(eval_diss_threshold)
}

getFoldDifferenceUpperBound<<-function()
{
  return(fold_diff_upper_bound)
}


getFoldDifferenceLowerBound<<-function()
{
  return(fold_diff_lower_bound)
}

getMinMotifLength<<-function()
{
  return(min_motif_length)
}



getVariableParameter<<-function()
{
  return(variable_parameter)
}

getDatasetDir<<-function()
{
  return(dataset_dir)
}

getRunsDir<<-function()
{
  return(runs_dir)
}

getOutputText<<-function()
{
  return(OutputText)
}

getOutputPlots<<-function()
{
  return(OutputPlots)
}

########################
actual_motifs_dir=file.path(runs_dir,"ActualMotifs")
dir.create(actual_motifs_dir, showWarnings = FALSE)

getActualMotifsDir<<-function()
{
  return(actual_motifs_dir)
}

########################
training_all_output_dir=file.path(runs_dir,"Training")
dir.create(training_all_output_dir, showWarnings = FALSE)

getTrainingAllOutputDir<<-function()
{
  return(training_all_output_dir)
}

########################
testing_all_output_dir=file.path(runs_dir,"Testing")
dir.create(testing_all_output_dir, showWarnings = FALSE)

getTestingAllOutputDir<<-function()
{
  return(testing_all_output_dir)
}

########################
runtime_evaluation_output_dir=file.path(runs_dir,"Evaluation of Runtime")
dir.create(runtime_evaluation_output_dir, showWarnings = FALSE)

getCompleteDataRuntimeEvaluationDir<<-function()
{
  return(runtime_evaluation_output_dir)
}


########################
performance_evaluation_output_dir=file.path(runs_dir,"Evaluation of Performance")
dir.create(performance_evaluation_output_dir, showWarnings = FALSE)

getTrainingTestingPerformanceEvaluationDir<<-function()
{
  return(performance_evaluation_output_dir)
}


########################
motifs_output_dir=file.path(runs_dir,"Predicted Motifs")
dir.create(motifs_output_dir, showWarnings = FALSE)

getMotifsOutputDir<<-function()
{
  return(motifs_output_dir)
}


########################
training_motifs_output_dir=file.path(motifs_output_dir,"Training")
dir.create(training_motifs_output_dir, showWarnings = FALSE)

getTrainingMotifsOutputDir<<-function()
{
  return(training_motifs_output_dir)
}

########################
testing_motifs_output_dir=file.path(motifs_output_dir,"Testing")
dir.create(testing_motifs_output_dir, showWarnings = FALSE)

getTestingMotifsOutputDir<<-function()
{
  return(testing_motifs_output_dir)
}

########################
unseen_testing_motifs_output_dir=file.path(motifs_output_dir,"UnseenTesting")
dir.create(unseen_testing_motifs_output_dir, showWarnings = FALSE)

getUnseenTestingMotifsOutputDir<<-function()
{
  return(unseen_testing_motifs_output_dir)
}


extractListParams<<-function(output_file_name_long)
{
  title_elements=unlist(strsplit(output_file_name_long,"_"))
  if(variable_parameter=="NumSeq")
  {
    title=paste(title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
                title_elements[length(title_elements)-3],"_",
                title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  } else if(variable_parameter=="SeqLength"){
    title=paste(title_elements[1],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
                title_elements[length(title_elements)-3],"_",
                title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
    
  } else if(variable_parameter=="NumClasses"){
    title=paste("C*300","_",title_elements[2],"\n",title_elements[length(title_elements)-4],"_",
                title_elements[length(title_elements)-3],"_",
                title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  }else if(variable_parameter=="MotifsPerClass"){
    title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",
                title_elements[length(title_elements)-3],"_",
                title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  } else if(variable_parameter=="MotifLength"){
    title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4]
                ,"_",
                title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  }else if(variable_parameter=="ProbaNoneOcc"){
    title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
                title_elements[length(title_elements)-3],"_",
                "_",title_elements[length(title_elements)],sep="")
  }else if(variable_parameter=="MutationRates"){
    title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
                title_elements[length(title_elements)-3],"_",
                title_elements[length(title_elements)-1],sep="")
  }
  return(title)
}

getParamsForPlotTitle<<-function()
{
  output_file_name_long=output_file_names_long[length(output_file_names_long)]
  title=extractListParams(output_file_name_long)
  return(title)
}







###############################################
# Generating Datasets
###############################################

all_datasets=list.files(datasets_dir,full.names=TRUE)#all input data files
output_file_names_short=c()
output_file_names_long=c()
seq_lengths_all_datasets=c()

for(i in 1:length(all_datasets))
{
  curr_dir=all_datasets[i]
  setwd(curr_dir)
  print(paste("WE ARE GENERATING DATASET ",i,"/",length(all_datasets),sep=""))  
  source(file.path(home_dir,"DataGenerator.R"))
  output_file_names_short=c(output_file_names_short,getOutputFileName())
  output_file_names_long=c(output_file_names_long,getOutputFileNameLong())
  seq_lengths_all_datasets=c(seq_lengths_all_datasets,getSequenceLength())
}

getOutputFileNamesShort<<-function()
{
  return(output_file_names_short)
}

getOutputFileNamesLong<<-function()
{
  return(output_file_names_long)
}

getSeqLengthsAllDatasets<<-function()
{
  return(seq_lengths_all_datasets)
}



###############################################
# Running the algorithm
###############################################

  data_type_vector=c()

  for(r in 1:num_runs)
  {
    print(paste("WE ARE AT RUN ",r,"/",num_runs,sep=""))
    
    current_run_name=paste("Run",r,sep="")
    getCurrentRunName<<-function()
    {
      return(current_run_name)
    }
    
    current_run_dir=file.path(runs_dir,paste("Run",r,sep=""))
    dir.create(current_run_dir, showWarnings = FALSE)
    
    #########################
    run_time_evaluation_dir=file.path(current_run_dir,"RunTimeEvaluation")
    dir.create(run_time_evaluation_dir, showWarnings = FALSE)
    
    getRunTimeEvaluationDir<<-function()
    {
      return(run_time_evaluation_dir)
    }
    
    getCurrentRunDir<<-function()
    {
      return(current_run_dir)
    }
    
    getAllOutputDir<<-function()
    {
      return(current_run_dir)
    }
    
          
    for(i in 1:length(all_datasets))
    {
      print(paste("WE ARE AT DATASET ",i,"/",length(all_datasets),sep=""))
      curr_dir=all_datasets[i]
      setwd(curr_dir)
      
      
      current_output_file_name_short<<-output_file_names_short[i]
      getOutputFileNameShort<<-function()
      {
        return(output_file_names_short[i])
      }
      
      getOutputFileNameLong<<-function()
      {
        return(output_file_names_long[i])
      }
      
      #importing current UNPARTITIONED datasets
      unpartitioned_dataset_filenames=c()
      unpartitioned_dataset_filenames <- list.files(file.path(curr_dir,"DataSets","with_motifs"), pattern="*.all", full.names=TRUE)
      unpartitioned_dataset_data_file_name=unpartitioned_dataset_filenames[1]
      unpartitioned_dataset_data = read.table(unpartitioned_dataset_data_file_name, quote="\"")
      print("dim(unpartitioned_dataset_data)")
      print(dim(unpartitioned_dataset_data))
      class_col=unpartitioned_dataset_data[,ncol(unpartitioned_dataset_data)]
      print("class_col")
      print(class_col)
      training_partition_indexes=unlist(createDataPartition(class_col,1,0.9))
      print("length(training_partition_indexes)")
      print(length(training_partition_indexes))
      training_partition=unpartitioned_dataset_data[training_partition_indexes,]
      print("dim(training_partition)")
      print(dim(training_partition))
      testing_partition=unpartitioned_dataset_data[-training_partition_indexes,]
      print("dim(testing_partition)")
      print(dim(testing_partition))

      #importing current testing datasets
#       unseen_testing_filenames=c()
#       unseen_testing_filenames <- list.files(file.path(curr_dir,"DataSets","with_motifs"), pattern="*.test", full.names=TRUE)
#       unseen_testing_data_file_name=unseen_testing_filenames[1]
#       unseen_testing_data = read.table(unseen_testing_data_file_name, quote="\"")


      #importing current dataset
#       filenames=c()
#       filenames <- list.files(file.path(curr_dir,"DataSets","with_motifs"), pattern="*.data", full.names=TRUE)
#       data_file_name=filenames[1]
#       data = read.table(data_file_name, quote="\"")

      data=training_partition
      getData<<-function()
      {
        return(data)
      }
      
      
      class_lables=data[,ncol(data)]
      data_type_vector=c()
      
      #create folds
      folds_indixes=createFolds(class_lables,num_folds)
#       
      for(f in 1:length(folds_indixes))
      {
        print(paste("WE ARE AT FOLD ",f,"/",length(folds_indixes),sep=""))
        
        current_fold_indixes=folds_indixes[[f]]
        print("length(current_fold_indixes)")
        print(length(current_fold_indixes))
        print("dim(data)")
        print(dim(data))
        data_type_vector=c("Training","Testing","UnseenTesting")
        current_fold_name=paste("Fold",f,sep="")
        
        
        for(t in 1:length(data_type_vector))
        {
          data_type<<-data_type_vector[t]
          
          getDataType<<-function()
          {
            return(data_type)
          }
          
          current_data=matrix()
          if(data_type=="Training")
          {
            current_data=data[-current_fold_indixes,]
#             data_type="Training"
          }
          else if(data_type=="Testing")
          {
            current_data=data[current_fold_indixes,]
#             data_type="Testing"
          }else if(data_type=="UnseenTesting"){
#             importing current UNSEEN TESTING dataset
                  current_data=testing_partition
          }


          print(paste("We are at ",data_type," dataset ",current_fold_name,sep=""))
#           current_data=unlist(complimentary_datasets)[t]
          print("dim(current_data)")
          print(dim(current_data))
          
          getCurrentFoldName<<-function()
          {
            return(current_fold_name)
            
          }
          
          
          getCurrentFold<<-function()
          {
            return(current_data)
          }
          
          
          #   ###########################
          current_dataset_type=file.path(current_run_dir,data_type)
          dir.create(current_dataset_type, showWarnings = FALSE)
          
          getCurrentDataSetTypeDir<<-function()
          {
            return(current_dataset_type)
          }
          
          #   ###########################
          current_fold_output_dir=file.path(current_dataset_type,current_fold_name)
          dir.create(current_fold_output_dir, showWarnings = FALSE)
          
          getCurrentFoldOutputDir<<-function()
          {
            return(current_fold_output_dir)
          }
          
          
          
          #   ###########################
          current_fold_dir=file.path(current_fold_output_dir,"1-Fold")
          dir.create(current_fold_dir, showWarnings = FALSE)
          
          getCurrentFoldsDir<<-function()
          {
            return(current_fold_dir)
          }
          

#           if(data_type!="UnseenTesting")
            write.table(current_data,paste(current_fold_dir,"/",current_fold_name,"_",current_output_file_name_short,".txt",sep=""))
#             write.table(current_fold_indixes,paste(current_fold_dir,"/",current_fold_name,"_",current_output_file_name_short,"_indixes.txt",sep=""))
    
          
          #########################
          allele_freq_scatter_plots_dir=file.path(current_fold_output_dir,"2-AlleleFreqScatterPlots")
          dir.create(allele_freq_scatter_plots_dir, showWarnings = FALSE)
          
          getAlleleFreqScatterPlotsDir<<-function()
          {
            return(allele_freq_scatter_plots_dir)
          }
          
          #####################
          pvalues_data_per_allele_dir=file.path(current_fold_output_dir,"3-PvaluesDataPerAllele")
          dir.create(pvalues_data_per_allele_dir, showWarnings = FALSE)
          
          getPvaluesDataPerAlleleDir<<-function()
          {
            return(pvalues_data_per_allele_dir)
          }
          
          #########################
          pvalue_scatter_plots_dir=file.path(current_fold_output_dir,"4-PvalueScatterPlots")
          dir.create(pvalue_scatter_plots_dir, showWarnings = FALSE)
          
          getPvalueScatterPlotsDir<<-function()
          {
            return(pvalue_scatter_plots_dir)
          }
          
          #########################
          pvalue_box_plots_dir=file.path(current_fold_output_dir,"5-PvalueBoxPlots")
          dir.create(pvalue_box_plots_dir, showWarnings = FALSE)
          
          getPvalueBoxPlotsDir<<-function()
          {
            return(pvalue_box_plots_dir)
          }
          
          #########################
          predicted_differentiable_areas_dir=file.path(current_fold_output_dir,"6-PredictedDifferentiableAreas")
          dir.create(predicted_differentiable_areas_dir, showWarnings = FALSE)
          
          getPredictedDifferentiableAreasDir<<-function()
          {
            return(predicted_differentiable_areas_dir)
          }
          
          ########################
          predicted_motifs_dir=file.path(current_fold_output_dir,"7-PredictedMotifs")
          dir.create(predicted_motifs_dir, showWarnings = FALSE)
          
          getPredictedMotifsDir<<-function()
          {
            return(predicted_motifs_dir)
          }
          
          ########################
          true_positive_motifs_dir=file.path(current_fold_output_dir,"8-MatchedMotifs")
          dir.create(true_positive_motifs_dir, showWarnings = FALSE)
          
          getTruePositiveMotifsDir<<-function()
          {
            return(true_positive_motifs_dir)
          }
          
          ########################
          trees_dir=file.path(current_fold_output_dir,"9-Trees")
          dir.create(trees_dir, showWarnings = FALSE)
          
          getTreesDir<<-function()
          {
            return(trees_dir)
          }
          
          #########################
          predicion_evaluation_dir=file.path(current_fold_output_dir,"10-PredictionEvaluation")
          dir.create(predicion_evaluation_dir, showWarnings = FALSE)
          
          getPredicionEvaluationDir<<-function()
          {
            return(predicion_evaluation_dir)
          }

          #########################
          run_time_evaluation_dir_per_fold=file.path(current_fold_output_dir,"11-RunTimeEvaluationPerFold")
          dir.create(run_time_evaluation_dir_per_fold, showWarnings = FALSE)
          getRunTimeEvaluationPerFoldDir<<-function()
          {
            return(run_time_evaluation_dir_per_fold)
          }
          
          #########################
          splsda_predicted_diff_locations_dir=file.path(current_fold_output_dir,"12-SplsdaPredictedDiffLocations")
          dir.create(splsda_predicted_diff_locations_dir, showWarnings = FALSE)
          
          getSplsdaPredictedDiffLocations<<-function()
          {
            return(splsda_predicted_diff_locations_dir)
          }
          
          if(OutputPlots)
          {
            source(file.path(home_dir,"PlotGenerator.R"))
          }
          source(file.path(home_dir,"RunSPLSDA.R"))
          source(file.path(home_dir,"AlgorithmSrc","RunAlgorithm.R"))
          if(OutputRuns)
          {
            source(file.path(home_dir,"RunSPLSDAforTiming.R"))
            source(file.path(home_dir,"AlgorithmSrc","RunAlgorithmForTiming.R"))
          }
        }#for(t in 1:length(data_type_vector))
      }#End of for(f in 1:num_folds)



      data_type<<-"complete"
#       source(file.path(home_dir,"AlgorithmSrc","RunAlgorithmForTiming.R"))
#       source(file.path(home_dir,"RunSPLSDAforTiming.R"))
    }#End of for(i in 1:num_datasets)

  }





##################################################
# Evaluation
##################################################

#########################
for(r in 1:num_runs)
{
  current_run_name=paste("Run",r,sep="")
  getCurrentRunName<<-function()
  {
    return(current_run_name)
  }
  
  current_run_dir=file.path(runs_dir,paste("Run",r,sep=""))
  dir.create(current_run_dir, showWarnings = FALSE)
    
    getAlleleFreqScatterPlotsDir<<-function(data_type_dir)
    {
    
      allele_freq_scatter_plots_dir=file.path(data_type_dir,"2-AlleleFreqScatterPlots")  
      return(allele_freq_scatter_plots_dir)
    }
    
    #####################
    
    getPvaluesDataPerAlleleDir<<-function(data_type_dir)
    {
      pvalues_data_per_allele_dir=file.path(data_type_dir,"3-PvaluesDataPerAllele")
      return(pvalues_data_per_allele_dir)
    }
    
    #########################
    
    getPvalueScatterPlotsDir<<-function(data_type_dir)
    {
      pvalue_scatter_plots_dir=file.path(data_type_dir,"4-PvalueScatterPlots")   
      return(pvalue_scatter_plots_dir)
    }
    
    #########################
    getPvalueBoxPlotsDir<<-function(data_type_dir)
    {
      pvalue_box_plots_dir=file.path(data_type_dir,"5-PvalueBoxPlots")
      return(pvalue_box_plots_dir)
    }
    
    #########################
    getPredictedDifferentiableAreasDir<<-function(data_type_dir)
    {
      predicted_differentiable_areas_dir=file.path(data_type_dir,"6-PredictedDifferentiableAreas")      
      return(predicted_differentiable_areas_dir)
    }
    
    ########################
    getPredictedMotifsDir<<-function(data_type_dir)
    {
      predicted_motifs_dir=file.path(data_type_dir,"7-PredictedMotifs")      
      return(predicted_motifs_dir)
    }
    
    ########################
    getTruePositiveMotifsDir<<-function(data_type_dir)
    {
      true_positive_motifs_dir=file.path(data_type_dir,"8-MatchedMotifs")      
      return(true_positive_motifs_dir)
    }
    
    ########################
    getTreesDir<<-function(data_type_dir)
    {
      trees_dir=file.path(data_type_dir,"9-Trees")
      return(trees_dir)
    }
    
    #########################
    getPredicionEvaluationDir<<-function(data_type_dir)
    {
      predicion_evaluation_dir=file.path(data_type_dir,"10-PredictionEvaluation")      
      return(predicion_evaluation_dir)
    }
    
    #########################
    getRunTimeEvaluationPerFoldDir<<-function(data_type_dir)
    {
      run_time_evaluation_dir_per_fold=file.path(data_type_dir,"11-RunTimeEvaluationPerFold")      
      return(run_time_evaluation_dir_per_fold)
    }
    
    #########################      
    getSplsdaPredictedDiffLocations<<-function(data_type_dir)
    {
      splsda_predicted_diff_locations_dir=file.path(data_type_dir,"12-SplsdaPredictedDiffLocations")
      return(splsda_predicted_diff_locations_dir)
    }
    

    for(f in 1:num_folds)
    {
      print(paste("EVALUATING PERFORMANCE. WE ARE AT FOLD ",f,"/",num_folds,sep=""))
      current_fold_name=paste("Fold",f,sep="")
      
      getCurrentFoldName<<-function()
      {
        return(current_fold_name)
        
      }
      
      training_dir=file.path(current_run_dir,"Training")
      current_fold_training_dir=file.path(training_dir,current_fold_name)
      
      getTrainingDir<<-function()
      {
        return(current_fold_training_dir)
      }
      
      testing_data_type=c("Testing","UnseenTesting")
      
      for(tt in 1:length(testing_data_type))
      {
        current_testing_data_type=testing_data_type[tt]
        print("current_testing_data_type")
        print(current_testing_data_type)
        testing_dir=file.path(current_run_dir,current_testing_data_type)
        current_fold_testing_dir=file.path(testing_dir,current_fold_name)
        
        getTestingDir<<-function()
        {
          return(current_fold_testing_dir)
        }
        
        getEvaluationOutputDir<<-function()
        {
          if(current_testing_data_type=="Testing")
            return(getTrainingAllOutputDir())
          else
            return(getTestingAllOutputDir())
        }        
        source(file.path(home_dir,"EvaluationSplsda.R"))
        source(file.path(home_dir,"EvaluationDiffAreas.R"))             
        source(file.path(home_dir,"EvaluationFullDiffMotifs.R"))
      }
      
      
    }#end of for every fold
    



    data_type<<-"complete"
#     source(file.path(home_dir,"EvaluationRuntime.R"))
#     source(file.path(home_dir,"EvaluationSplsdaRuntime.R"))
}#End of for(r in 1:num_runs)





##################################################
# Generation of plots
##################################################


data_type_vector<<-c("Training","Testing")
for(t in 1:length(data_type_vector))
{
  data_type<<-data_type_vector[t]
  source(file.path(home_dir,"SummarySplsdaEvaluation.R"))
}

source(file.path(home_dir,"SummaryDiffAreasEvaluation.R"))
source(file.path(home_dir,"SummaryFullDiffMotifsEvaluation.R"))

# data_type_vector=c("Testing")
for(t in 1:length(data_type_vector))
{
  data_type<<-data_type_vector[t]
  source(file.path(home_dir,"PlottingRocCurves.R"))
}

data_type<<-"complete"
# source(file.path(home_dir,"SummaryRuntimeEvaluation.R"))
# source(file.path(home_dir,"SummarySplsdaRuntimeEvaluation.R"))




time=proc.time() - ptm  
print("ALL DONE!!!")
print(paste("TIME to run all runs on all datasets: ",time['elapsed']))
