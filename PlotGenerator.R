#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Plot Generator: 
# A program that generates plots and tables to analyse randomly generated genome data
# Author: Gioia Wehbe
# LAST UPDATED: 31-01-2015
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print("Inside PlotGenerator.R")

# Start the clock!
ptm <- proc.time()
date="31-01-2015"

actual_motifs_locations_dir=getActualMotifsDir()
########################################################################################
#loading data
########################################################################################
print("Plotting p-values figures")
# #getting dataset file name

current_fold=getCurrentFold()

OutputPlots=getOutputPlots()

#Loading Parameters table
parameter_table = read.table(paste(getwd(),"/Parameters/parameters.txt",sep=""),header=FALSE,sep="\t",strip.white=TRUE)
parameter_table=as.matrix(parameter_table)  #import input file containing parameters as a table and convert it to a matrix
num_classes=as.numeric(parameter_table[1,2])            #number of classes specified by the user in the parameter file
num_SNPs= as.numeric(parameter_table[3,2])             #sequence length specified by the user in the parameter file
class_probas_text=parameter_table[4,2]  #A string containing list of class probas
class_probas_list=strsplit(class_probas_text, ",") #list of class probas
class_probas_vector=class_probas_list[[1]]        # convert list to vector
nuc_probas_text=parameter_table[5,2]  #A string containing list of nucleotide probas
nuc_probas_list=strsplit(nuc_probas_text, ",")     # stores neucleotide probas as a list
nuc_probas_vector=nuc_probas_list[[1]]        # convert list to vector
output_file_name=parameter_table[6,2]   #output files name base on variable parameters
output_file_name_long=parameter_table[7,2]   #long output files name base on input parameters

nuc_list=c('0','1','2')                   #list of genome nucleotides
colnames(current_fold)=append(paste(rep("SNP",num_SNPs),1:num_SNPs),"Class")

#Import Locations Table and then fix column names (colClasses="character" is the trick):

actual_motifs_file_name=paste(actual_motifs_locations_dir,"/ActualMotifLocations_",output_file_name,".csv",sep="")
current_fold_Motif_Locations <- read.csv(actual_motifs_file_name,header=TRUE,sep=",",colClasses="character")
rownames(data_Motif_Locations)=NULL
data_Motif_Locations=t(data_Motif_Locations)
data_Motif_Locations=data_Motif_Locations[-1,]
num_motifs=ncol(data_Motif_Locations)

#functions that return the list of start/end locations given a certain class

getMotifsStartLocations<-function(class)
{
  start_locations_list=c()
  for(i in 1:ncol(data_Motif_Locations))
  {
    if(data_Motif_Locations[1,i]==class)
    {
      start_locations_list= append(start_locations_list,data_Motif_Locations[3,i])
    }
  }
  return(start_locations_list)
}

getMotifsEndLocations<-function(class)
{
  end_locations_list=c()
  for(i in 1:ncol(data_Motif_Locations))
  {
    if(data_Motif_Locations[1,i]==class)
    {
      end_locations_list= append(end_locations_list,data_Motif_Locations[4,i])
    }
  }
  return(end_locations_list)
}

#Get motif probas of nonoccurence and motif mutation rates in lists

motif_start_locations=as.numeric(data_Motif_Locations[3,1:ncol(data_Motif_Locations)])
row.names(motif_start_locations)=NULL
classes_per_motif=as.numeric(data_Motif_Locations[1,1:ncol(data_Motif_Locations)])
row.names(classes_per_motif)=NULL

#Sort the data set by the class label and put each class in a separate data set
data_per_class_list=list(1:num_classes)#NB!! it should be a list
data_per_class_size_list=c(1:num_classes)
for(i in 1:num_classes)
{
  class_index=which(current_fold$Class==i)
  class_seqs=current_fold[class_index,]
  data_per_class=as.data.frame(class_seqs)
  data_per_class_size=nrow(data_per_class)
  data_per_class$Class=NULL
  rownames(data_per_class)=NULL
  data_per_class_list[[i]]=data_per_class
  data_per_class_size_list[[i]]=data_per_class_size
}

########################################################################################
#Functions used to create frequency tables
########################################################################################

generateFreqTableRow<-function(class,SNP_type,class_data,num_SNPs)
{
  freqs=rep(NA,num_SNPs+2)
  freqs[1]=class
  freqs[2]=SNP_type
  for(i in 1:num_SNPs)
  {
    snp=class_data[,i]
    snp_table=table(snp)    
    all_freq=as.data.frame(snp_table)
    rownames(all_freq)=all_freq[,1]
    
    SNP_type_ind=as.character(SNP_type)
    SNP_type_freq=all_freq[SNP_type_ind,2]
    if(is.na(SNP_type_freq))
      SNP_type_freq=0
    freqs[i+2]=SNP_type_freq
  }
  return(freqs)
}

generateFreqTable<-function(num_classes,num_SNPs,class_data_list)#class data list is a list of sequences per class
{
  
  freq_table=as.data.frame(matrix(nrow=num_classes*3,ncol=num_SNPs+2))
  row_index=1
  for(c in 1:num_classes)
  {
    class_data=as.data.frame(class_data_list[c])
    for(s in 0:2)
    {
      row=generateFreqTableRow(c,s,class_data,num_SNPs)
      freq_table[row_index,]=row
      row_index=row_index+1
    }
  }
  return(freq_table)
}

#############################
# Generate Frequency table
#############################

data_freq_table=generateFreqTable(num_classes,num_SNPs,data_per_class_list)
colnames(data_freq_table)=append(c("Class","Allele"),colnames(data[1:num_SNPs]))

####################################
# Generate Relative Frequency table
####################################

data_rel_freq_table=data.frame()
class_first_row=1
for(i in 1:num_classes)
{
  data_rel_freq_table=rbind(data_rel_freq_table,data_freq_table[class_first_row,]/data_per_class_size_list[i]
                            ,data_freq_table[class_first_row+1,]/data_per_class_size_list[i]
                            ,data_freq_table[class_first_row+2,]/data_per_class_size_list[i])
  class_first_row=class_first_row+3
}

data_rel_freq_table[,1]=data_freq_table[,1]
data_rel_freq_table[,2]=data_freq_table[,2]
colnames(data_rel_freq_table)=colnames(data_freq_table)


#####################################################
# Initializing ploting constants
#####################################################

class_color_vector=c("red","blue" ,"yellow", "green"  ,"purple" ,"grey" ,"pink","cyan","magenta")

motif_color_vector=c("red4","blue4" ,"yellow4", "green4" ,"purple4" ,"grey4" ,"pink4","cyan4","magenta4")


#########################################################
# Extracting parameters for plot title
#########################################################

extractListParams<-function(output_file_name_long)
{
  title_elements=unlist(strsplit(output_file_name_long,"_"))
  title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
              title_elements[length(title_elements)-3],"_",title_elements[length(title_elements)-2],"_",
              title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  return(title)
}

getParamsForPlotTitle<-function()
{
  title=extractListParams(output_file_name_long)
  return(title)
}

#####################################################################################
#Plotting the relative frquencies of each allele in a certain class (Scatter Plot)
#####################################################################################


plotRelativeFreqOfAllelePerClass<-function(num_SNPs,class,allele,row,col_class,col_motif,motif_start_list,motif_end_list,
                                           allele_pch)
{
  
  par(mar=c(5.1,8,4,8),xpd=TRUE)
  
  plot(x=NULL, y=NULL, xlab="Position on chromosome", ylab=paste("Rel. Freq. of Allele",allele,"in class",class),
       xlim=c(1,num_SNPs),ylim=c(0,1))
  
  num_bins=num_classes
  colcode=class_color_vector[1:num_classes]
  legend("topright", inset=c(-0.4,0),legend = class_probas_vector,title = "Class Probabilities",fill = colcode,
         bty = "n") # border
  
  num_bins=length(nuc_probas_vector)
  colcode=class_color_vector[1:num_motifs]
  
  legend("right",inset=c(-0.4,0), legend = nuc_probas_vector,title = "Allele Probabilities",pch=nuc_list,
         bty = "n") # border
  
  title=getParamsForPlotTitle()
  title(main=title)
  
  rel_freqs=data_rel_freq_table[row,3:(num_SNPs+2)]
  points(x=1:num_SNPs, y=rel_freqs, cex=0.8, col=col_class,pch=allele_pch)
  for(i in 1:length(motif_start_list))
  {
    motif_start=motif_start_list[i]
    motif_end=motif_end_list[i]
    points(x=motif_start:motif_end, y=rel_freqs[motif_start:motif_end], cex=0.8, col=col_motif,pch=allele_pch)
  }
  
}




if(OutputPlots)
{
  fname=paste(getAlleleFreqScatterPlotsDir(),"/RelFreqPlots_",output_file_name,".pdf",sep="")
  pdf(fname)
  row=1
  for(c in 1:num_classes)
  {
    col_class=class_color_vector[c]
    col_motif=motif_color_vector[c]
    motif_start_list=getMotifsStartLocations(c)
    motif_end_list=getMotifsEndLocations(c)
    for(a in 0:2)
    {
      allele_pch=as.character(a)
      plotRelativeFreqOfAllelePerClass(num_SNPs,c,a,row,col_class,col_motif,motif_start_list,motif_end_list,allele_pch)
      row=row+1
    }
  }
  dev.off ();
}
  


####################################################################################
# compute and plot the p-values of every allele for all the classes
####################################################################################
# 

computeSnpPvalue<-function(snp,allele)
{  
  class_matrix=data.frame(row.names=NULL)
  row=allele+1
  for(c in 1:num_classes)#TO DO: Fix by using which...
  {
    class_successes=as.numeric(data_freq_table[row,snp+2])
    class_failures=as.numeric(data_per_class_size_list[[c]]-class_successes)
    class_matrix=rbind(class_matrix,list(class_successes,class_failures))
    row=row+3
  }
  rownames(class_matrix)=NULL
  p_value=prop.test(as.matrix(class_matrix))$p.value
  p_value=log10(p_value)
  p_value=abs(p_value)
  return(p_value)
}

plotPvalues<-function(p_values,allele,num_SNPs,num_classes,class_color_vector)
{
  
  par(mar=c(5.1,8,4,8),xpd=TRUE)
  plot(x=NULL, y=NULL, xlab="Position on chromosome", ylab=paste("abs(log(p-values)) of allele",allele),xlim=c(1,num_SNPs),ylim=c(0,500))
  
  num_bins=num_classes
  colcode=class_color_vector[1:num_classes]
  legend("topright", inset=c(-0.4,0),legend = class_probas_vector,title = "Class Probabilities",fill = colcode,
         bty = "n") # border
  
  #allele legend
  num_bins=length(nuc_probas_vector)
  colcode=class_color_vector[1:num_motifs]
  
  legend("right",inset=c(-0.4,0), legend = nuc_probas_vector,title = "Allele Probabilities",pch=nuc_list,
         bty = "n") # border
title_elements=unlist(strsplit(output_file_name_long,"_"))
  title=paste(title_elements[1],"_",title_elements[2],"_",title_elements[3],"\n",title_elements[length(title_elements)-4],"_",
              title_elements[length(title_elements)-3],"_",title_elements[length(title_elements)-2],"_",
              title_elements[length(title_elements)-1],"_",title_elements[length(title_elements)],sep="")
  title(main=title)
  
  points(x=1:num_SNPs, y=p_values, cex=0.8, col="wheat3" ,pch=as.character(allele))
  for(c in 1:num_classes)
  {
    class=c
    motif_start_list=getMotifsStartLocations(class)
    motif_end_list=getMotifsEndLocations(class)
    for(i in 1:length(motif_start_list))
    {
      col_motif=class_color_vector[as.numeric(class)]
      motif_start=motif_start_list[i]
      motif_end=motif_end_list[i]
      points(x=motif_start:motif_end, y=p_values[motif_start:motif_end], cex=0.8, col=col_motif,pch=as.character(allele))
    }
  }
}

if(OutputPlots)
{
  pdf(paste(getPvalueScatterPlotsDir(),"/PvaluePlots_",output_file_name,".pdf",sep="")) 
}

file_pvalues_allel0=paste(getPvaluesDataPerAlleleDir(),"/PvaluesAllele0.txt",sep="")
file_pvalues_allel1=paste(getPvaluesDataPerAlleleDir(),"/PvaluesAllele1.txt",sep="")
file_pvalues_allel2=paste(getPvaluesDataPerAlleleDir(),"/PvaluesAllele2.txt",sep="")

p_values_matrix_all_alleles=matrix(nrow=1,ncol=num_SNPs+1)
for(a in 0:2)
{
  p_values=mapply(computeSnpPvalue,1:num_SNPs,MoreArgs=list(a))
  p_values_row=append(output_file_name,p_values)
  p_values_matrix=matrix(p_values_row,nrow=1,byrow=TRUE)
  p_values_matrix_all_alleles=rbind(p_values_matrix_all_alleles,append(a,p_values))
  if(a==0)
    write.table(p_values_matrix,file=file_pvalues_allel0,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  else if(a==1)
    write.table(p_values_matrix,file=file_pvalues_allel1,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  else if(a==2)
    write.table(p_values_matrix,file=file_pvalues_allel2,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  
  if(OutputPlots)
  {
    plotPvalues(p_values,a,num_SNPs,num_classes,class_color_vector)
  }
}
if(OutputPlots)
{
  dev.off (); 
}



p_values_matrix_all_alleles=p_values_matrix_all_alleles[-1,]
rownames(p_values_matrix_all_alleles)=NULL
cnames=colnames(data_freq_table)
cnames=cnames[2:length(cnames)]
colnames(p_values_matrix_all_alleles)=cnames


getPvaluesMatrixAllAlleles<<-function()
{
  return(p_values_matrix_all_alleles)
}


###################
# Output runtime
###################

time=proc.time() - ptm            #compute the time taken by the script

print("PLOT GENERATING CODE ENDED! Plots have been saved in output folders.")   #print out in terminal an indication that the code has ended 
print("Runtime (in seconds)")
print(time['elapsed']) # Stop the clock and print out time elapse in the terminal 
