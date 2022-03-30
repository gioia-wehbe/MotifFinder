#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# InputGenerator: 
# A program that generates input files and directories to be used by the DataGenerator taking into consideration all possible 
# combinations of the parameter possibilities specified by the user. To generate the set of datasets, put the 
# InputGenerator.R source file in the folder where you want the datasets to be generated and then run the R script from cmd
# passing the number of runs as a parameter
# Author: Gioia Wehbe
# LAST UPDATED: 10-06-2015 (Fixed names and passes variable parameter and path as cmd arguments)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

args <- commandArgs(trailingOnly = TRUE)
run_num=as.numeric(args[1])
vriable_parameter=args[2]
output_path=args[3]

#for testing:
# run_num=1
# vriable_parameter="MotifLength"


allele_probas_poss=list(c(0.8,0.15,0.05))
motif_num_occ_per_seq_poss=c(1)

#Best Case:
num_seq_poss=c(6000)    #FIXED
seq_length_poss=c(8000)   #FIXED
num_classes_poss=c(6)   #FIXED
class_probas_poss_2=list(c(0.8,0.2))
class_probas_poss_3=list(c(0.5,0.3,0.2))
class_probas_poss_5=list(c(0.15,0.15,0.20,0.30,0.20))
class_probas_poss_6=list(c(0.15,0.15,0.20,0.20,0.20,0.10))
class_probas_poss_8=list(c(0.1,0.1,0.2,0.2,0.1,0.1,0.1,0.1))
class_probas_poss_10=list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
class_probas_poss_14=list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05))
max_motifs_per_class_poss=c(1)  #FIXED
motif_lengths_poss=c(5)  #FIXED
motif_proba_of_none_occ_poss=c(0.1)    #FIXED
motif_mutation_rates_poss=c(0.03) #FIXED


#Medians
# num_seq_poss=c(3000)    #FIXED
# seq_length_poss=c(6000)   #FIXED
# num_classes_poss=c(3)   #FIXED
# class_probas_poss_2=list(c(0.8,0.2))
# class_probas_poss_3=list(c(0.5,0.3,0.2))
# class_probas_poss_5=list(c(0.15,0.15,0.20,0.30,0.20))
# class_probas_poss_6=list(c(0.15,0.15,0.20,0.20,0.20,0.10))
# class_probas_poss_8=list(c(0.1,0.1,0.2,0.2,0.1,0.1,0.1,0.1))
# class_probas_poss_10=list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
# class_probas_poss_14=list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05))
# max_motifs_per_class_poss=c(3)  #FIXED
# motif_lengths_poss=c(15)  #FIXED
# motif_proba_of_none_occ_poss=c(0.15)    #FIXED
# motif_mutation_rates_poss=c(0.015) #FIXED



if(vriable_parameter=="NumSeq") 
  num_seq_poss=c(1000,2000,3000,4000,5000) else
if(vriable_parameter=="SeqLength")
  seq_length_poss=c(2000,4000,8000,16000) else
if(vriable_parameter=="NumClasses")
{
  num_classes_poss=c(2,6,10,14) #don't forget to fix if else statements if numclasses possibilties change
  class_probas_poss_2=list(c(0.5,0.5))
  p=1/3
  class_probas_poss_3=list(c(p,p,p))
  p=1/5
  class_probas_poss_5=list(c(p,p,p,p,p))
  p=1/6
  class_probas_poss_6=list(c(p,p,p,p,p,p))
  p=1/8
  class_probas_poss_8=list(c(p,p,p,p,p,p,p,p))
  p=1/10
  class_probas_poss_10=list(c(p,p,p,p,p,p,p,p,p,p))
  p=1/14
  class_probas_poss_14=list(c(p,p,p,p,p,p,p,p,p,p,p,p,p,p))
} else
if(vriable_parameter=="MotifsPerClass")
  max_motifs_per_class_poss=c(1,2,3,4,5) else
if(vriable_parameter=="MotifLength")
{
  motif_lengths_poss=c(5,10,15,20,25)
}else
if(vriable_parameter=="ProbaNoneOcc")
{
  motif_proba_of_none_occ_poss=c(0.05,0.1,0.15,0.2,0.25,0.3)
}else
if(vriable_parameter=="MutationRates")
{
  motif_mutation_rates_poss=c(0.01,0.02,0.03,0.04,0.05)
}

# print(paste("motif_lengths_poss",motif_lengths_poss))

main_dir=output_path
data_set_dir<<-"1-SimulatedDataSets"
sample_dir<<-"Sample"
input_dir<<-"Parameters"
output_dir<<-"DataSets"
no_motifs_dir<<-"no_motifs"
with_motifs_dir<<-"with_motifs"
sample_counter=1

setupNewSampleDirs<-function(output_files_name_short,run_dir)
{
  sample_dir<<-file.path(run_dir,paste(sample_dir,"_",output_files_name_short,sep=""))
  dir.create(sample_dir, showWarnings = FALSE)
#   print("sample_dir")
#   print(sample_dir)
  input_dir<<-file.path(sample_dir,input_dir)
  dir.create(input_dir, showWarnings = FALSE)
#   print("input_dir")
#   print(input_dir)
  output_dir<<-file.path(sample_dir,output_dir)
  dir.create(output_dir, showWarnings = FALSE)
#   print("output_dir")
#   print(output_dir)
  no_motifs_dir<<-file.path(output_dir,no_motifs_dir)
  dir.create(no_motifs_dir, showWarnings = FALSE)
#   print("no_motifs_dir")
#   print(no_motifs_dir)
  with_motifs_dir<<-file.path(output_dir,with_motifs_dir)
  dir.create(with_motifs_dir, showWarnings = FALSE)
#   print("with_motifs_dir")
#   print(with_motifs_dir)
  
}

resetDirectories<-function()
{
  main_dir=output_path
#   data_set_dir="1-SimulatedDataSets"
  sample_dir<<-"Sample"
  input_dir<<-"Parameters"
  output_dir<<-"DataSets"
  no_motifs_dir<<-"no_motifs"
  with_motifs_dir<<-"with_motifs"
}

generateParametersInputFile<-function(curr_num_classes,curr_num_seq,curr_seq_length
                    ,curr_class_probas,curr_allele_probas,output_files_name_short,output_files_name)
{
  parameter_file=file.path(input_dir,"parameters.txt")
  write(paste("Number of classes\t",curr_num_classes,
              "\nNumber of sequences\t",curr_num_seq,
              "\nLength of sequences\t",curr_seq_length,
              "\nProbability of each class\t",paste(curr_class_probas,collapse=","),
              "\nProbability of each nucleotide (0,1,2)\t",paste(curr_allele_probas,collapse=","),
              "\nOutput file name short (without extension)\t",output_files_name_short,
              "\nOutput file name long (without extension)\t",output_files_name),
              parameter_file)
                          
}

generateNumMotifsPerClassVector<-function(curr_num_classes,curr_max_motifs_per_class)
{
  num_motifs_vec=rep(NA,curr_num_classes)
  for(i in 1:curr_num_classes)
  {
    num_motifs=curr_max_motifs_per_class
    num_motifs_vec[i]=num_motifs
  }
  return(num_motifs_vec)
}

# generateNumMotifsPerClassVector<-function(curr_num_classes,curr_max_motifs_per_class)
# {
#   num_motifs_vec=rep(NA,curr_num_classes)
#   for(i in 1:curr_num_classes)
#   {
#     num_motifs=sample(curr_max_motifs_per_class,size=1)
#     num_motifs_vec[i]=num_motifs
#   }
#   return(num_motifs_vec)
# }


motifAlreadyExists<-function(all_motifs,current_motif)
{
  for(i in 1:length(all_motifs))
  {
    curr_motif_list=all_motifs[i]
    if(is.na(curr_motif_list)==FALSE)
    {
      curr_motif_vec=unlist(strsplit(curr_motif_list,","))
      if(current_motif %in% curr_motif_vec)
        return(TRUE)
    }
    else
      return(FALSE)
  }
  return(FALSE)
}

generateMotifsCol<-function(curr_num_classes,curr_max_motifs_per_class,curr_motif_lengths,num_motifs_per_class_vec)
{
  motifs_col=rep(NA,curr_num_classes)
  for(i in 1:curr_num_classes)
  {
    num_motifs=num_motifs_per_class_vec[i]
    class_motifs=rep(NA,num_motifs)
    motif=""
    for(j in 1:num_motifs)
    {
      repeat
      {
        motif_vector=sample(c('0','1','2'),size=as.numeric(curr_motif_lengths),replace=TRUE)
        motif=paste(motif_vector, collapse="-")
        if(motifAlreadyExists(motifs_col,motif)==FALSE)
          break
        else
          print("Motif already exists!!!")
      }
      class_motifs[j]=motif
    }
    class_motifs_comma_delimitted=paste(class_motifs,collapse=",")
    motifs_col[i]=class_motifs_comma_delimitted
  }
  return(motifs_col)
}

generateNumOccCol<-function(curr_num_classes,curr_motif_num_occ_per_seq,num_motifs_per_class_vec)
{
  num_occ_vec=rep(NA,curr_num_classes)
  for(i in 1:curr_num_classes)
  {
    num_occ_per_motif_vect=rep(curr_motif_num_occ_per_seq,num_motifs_per_class_vec[i])
    num_occ_per_motif_comma_dellimitted=paste(num_occ_per_motif_vect,collapse=",")
    num_occ_vec[i]=num_occ_per_motif_comma_dellimitted
  }
  return(num_occ_vec)
}


generateProbNoneOccCol<-function(curr_num_classes,curr_motif_proba_of_none_occ,num_motifs_per_class_vec)
{
  proba_none_occ_vec=rep(NA,curr_num_classes)
  for(i in 1:curr_num_classes)
  {
    proba_none_occ_per_motif_vect=rep(curr_motif_proba_of_none_occ,num_motifs_per_class_vec[i])
    proba_none_occ_per_motif_comma_dellimitted=paste(proba_none_occ_per_motif_vect,collapse=",")
    proba_none_occ_vec[i]=proba_none_occ_per_motif_comma_dellimitted
  }
  return(proba_none_occ_vec)
}


generateMutationRatesCol<-function(curr_num_classes,curr_motif_mutation_rates,num_motifs_per_class_vec)
{
  mutation_rates_vec=rep(NA,curr_num_classes)
  for(i in 1:curr_num_classes)
  {
    mutation_rates_per_motif_vect=rep(curr_motif_mutation_rates,num_motifs_per_class_vec[i])
    mutation_rates_per_motif_comma_dellimitted=paste(mutation_rates_per_motif_vect,collapse=",")
    mutation_rates_vec[i]=mutation_rates_per_motif_comma_dellimitted
  }
  return(mutation_rates_vec)
}



generateMotifsInputFile<-function(curr_num_classes,curr_max_motifs_per_class,
          curr_motif_lengths,curr_motif_num_occ_per_seq,
          curr_motif_proba_of_none_occ,curr_motif_mutation_rates)
{
  motifs_file=file.path(input_dir,"motifs.csv")
  classes_col=c(1:curr_num_classes)
  num_motifs_per_class_vec=generateNumMotifsPerClassVector(curr_num_classes,curr_max_motifs_per_class)
  motifs_col=generateMotifsCol(curr_num_classes,curr_max_motifs_per_class,curr_motif_lengths,num_motifs_per_class_vec)
  num_occ_col=generateNumOccCol(curr_num_classes,curr_motif_num_occ_per_seq,num_motifs_per_class_vec)
  prob_none_occ_col=generateProbNoneOccCol(curr_num_classes,curr_motif_proba_of_none_occ,num_motifs_per_class_vec)
  mutation_rates_col=generateMutationRatesCol(curr_num_classes,curr_motif_mutation_rates,num_motifs_per_class_vec) 
  motifs_table=matrix(c(classes_col,motifs_col,num_occ_col,prob_none_occ_col,mutation_rates_col),nrow=curr_num_classes,ncol=5)
  
  
  colnames(motifs_table)=c("Class","List_of_Motifs","Max_Number_of_Occurrences_in_a_Sequence",
                           "Probability_None_Occurrence_in_a_Sequence","Mutation_Rate_Per_Position")
  write.csv(motifs_table,motifs_file,quote=TRUE,row.names=FALSE)
}







for(curr_run_num in 1:run_num)
{
  curr_run_dir=paste("Run",curr_run_num,sep="")
  run_dir<<-file.path(main_dir,curr_run_dir)
  dir.create(run_dir, showWarnings = FALSE)
#   print("run_dir")
#   print(run_dir)
  current_data_set_dir<<-file.path(run_dir,data_set_dir)
  dir.create(current_data_set_dir, showWarnings = FALSE)
#   print("current_data_set_dir")
#   print(current_data_set_dir)
  for(num_seq_index in 1:length(num_seq_poss))
  {
    curr_num_seq=num_seq_poss[num_seq_index]
    if(nchar(curr_num_seq)==3)
      curr_num_seq=paste("000",curr_num_seq,sep="")
    if(nchar(curr_num_seq)==4)
      curr_num_seq=paste("00",curr_num_seq,sep="")
    if(nchar(curr_num_seq)==5)
      curr_num_seq=paste("0",curr_num_seq,sep="")
    for(seq_length_index in 1:length(seq_length_poss))
    {
      curr_seq_length=seq_length_poss[seq_length_index]
      if(nchar(curr_seq_length)==3)
        curr_seq_length=paste("000",curr_seq_length,sep="")
      if(nchar(curr_seq_length)==4)
        curr_seq_length=paste("00",curr_seq_length,sep="")
      if(nchar(curr_seq_length)==5)
        curr_seq_length=paste("0",curr_seq_length,sep="")
      for(num_classes_index in 1:length(num_classes_poss))
      {
        curr_num_classes=num_classes_poss[num_classes_index]
        if(curr_num_classes==2)
          class_probas_poss=class_probas_poss_2
        else if(curr_num_classes==3)
          class_probas_poss=class_probas_poss_3
        else if(curr_num_classes==6)
          class_probas_poss=class_probas_poss_6
        else if(curr_num_classes==10)
          class_probas_poss=class_probas_poss_10
        else if(curr_num_classes==14)
          class_probas_poss=class_probas_poss_14
        if(length(num_classes_poss)>1)
        {
          curr_num_seq=curr_num_classes*300
          if(nchar(curr_num_seq)==3)
            curr_num_seq=paste("000",curr_num_seq,sep="")
          if(nchar(curr_num_seq)==4)
            curr_num_seq=paste("00",curr_num_seq,sep="")
          if(nchar(curr_num_seq)==5)
            curr_num_seq=paste("0",curr_num_seq,sep="")
        }
        for(class_probas_index in 1:length(class_probas_poss))
        {
          curr_class_probas=class_probas_poss[[class_probas_index]]
          for(allele_probas_index in 1:length(allele_probas_poss))
          {
            curr_allele_probas=allele_probas_poss[[allele_probas_index]]
            for(max_motifs_per_class_index in 1:length(max_motifs_per_class_poss))
            {
              curr_max_motifs_per_class=max_motifs_per_class_poss[max_motifs_per_class_index]
              for(motif_lengths_index in 1:length(motif_lengths_poss))
              {
                curr_motif_lengths=motif_lengths_poss[motif_lengths_index]
                for(motif_num_occ_per_seq_index in 1:length(motif_num_occ_per_seq_poss))
                {
                  curr_motif_num_occ_per_seq=motif_num_occ_per_seq_poss[motif_num_occ_per_seq_index]
                  for(motif_proba_of_none_occ_index in 1:length(motif_proba_of_none_occ_poss))
                  {
                    curr_motif_proba_of_none_occ=motif_proba_of_none_occ_poss[motif_proba_of_none_occ_index]
                    for(motif_mutation_rates_index in 1:length(motif_mutation_rates_poss))
                    {
                      curr_motif_mutation_rates=motif_mutation_rates_poss[motif_mutation_rates_index]
                      if(nchar(curr_motif_proba_of_none_occ)==3)
                        curr_motif_proba_of_none_occ=paste(curr_motif_proba_of_none_occ,"0",sep="")
                      if(nchar(curr_motif_mutation_rates)==4)
                        curr_motif_mutation_rates=paste(curr_motif_mutation_rates,"0",sep="")
                      if(nchar(curr_motif_lengths)==1)
                        curr_motif_lengths=paste("0",curr_motif_lengths,sep="")
                      if(nchar(curr_num_classes)==1)
                        curr_num_classes=paste("0",curr_num_classes,sep="")
                      output_files_name=paste("S",curr_num_seq,"_L",curr_seq_length,"_C",curr_num_classes,"_CP",
                                              round(curr_class_probas[1],digits=3),"_N",paste(curr_allele_probas,collapse="_"),
                                              "_MperC",curr_max_motifs_per_class,"_LM",curr_motif_lengths,"_MnumOcc",
                                              curr_motif_num_occ_per_seq,"_PnoneOcc",curr_motif_proba_of_none_occ,"_MutR",
                                              curr_motif_mutation_rates,sep="")
                      if(vriable_parameter=="NumSeq") 
                        output_files_name_short=paste("NumSeq_",curr_num_seq,sep="") else
                      if(vriable_parameter=="SeqLength")
                        output_files_name_short=paste("SeqLength_",curr_seq_length,sep="") else
                      if(vriable_parameter=="NumClasses")
                        output_files_name_short=paste("NumClasses_",curr_num_classes,sep="") else
                      if(vriable_parameter=="MotifsPerClass") 
                        output_files_name_short=paste("MotifsPerClass_",curr_max_motifs_per_class,sep="") else
                      if(vriable_parameter=="MotifLength")
                        output_files_name_short=paste("MotifLength_",curr_motif_lengths,sep="") else
                      if(vriable_parameter=="ProbaNoneOcc")
                        output_files_name_short=paste("ProbaNoneOcc_",curr_motif_proba_of_none_occ,sep="") else
                      if(vriable_parameter=="MutationRates")
                        output_files_name_short=paste("MutationRates_",curr_motif_mutation_rates,sep="")
                          
                      
#                       print(paste("output_files_name_short",output_files_name_short))
                      
#                       print("Before setupNewSampleDirs")
                      setupNewSampleDirs(output_files_name_short,current_data_set_dir)
#                       print("After setupNewSampleDirs")
                      sample_counter=sample_counter+1
                      generateParametersInputFile(curr_num_classes,curr_num_seq,curr_seq_length,curr_class_probas,curr_allele_probas,output_files_name_short,output_files_name)
                      generateMotifsInputFile(curr_num_classes,curr_max_motifs_per_class,
                                              curr_motif_lengths,curr_motif_num_occ_per_seq,
                                              curr_motif_proba_of_none_occ,curr_motif_mutation_rates)
                      resetDirectories()
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}











