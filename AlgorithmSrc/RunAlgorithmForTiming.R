###########################################
# RunAlgorithm
# Last updated: 23-02-15
###########################################

print("Inside RunAlgorithmForTiming.R")

home_dir=getHomeDir()
OutputPlots=getOutputPlots()


filenames=c()
filenames <- list.files(paste(getwd(),"/DataSets/with_motifs",sep=""), pattern="*.data", full.names=TRUE)
file_name=filenames[1]
# current_output_file_name_short=getOutputFileNameShort()
print("current_output_file_name_short")
print(current_output_file_name_short)
#############################################################
#Import data and then fix column names:
#############################################################

# data <- read.table(file_name, quote="\"")
current_fold=getCurrentFold()
if(data_type=="complete")
{
  print("inside if(data_type==complete)")
  current_fold=getData() #evaluate time on complete data set
  print("dim(current_fold)")
  print(dim(current_fold))
  
}
class_column=current_fold[,ncol(current_fold)]
current_fold_no_class_column=current_fold[,-length(current_fold)]
colnames(current_fold_no_class_column)=1:ncol(current_fold_no_class_column)
num_classes=length(levels(as.factor(current_fold[,ncol(current_fold)])))
num_sequences= nrow(current_fold)             
num_SNPs= ncol(current_fold)-1          
nuc_list=c('0','1','2')
colnames(current_fold)=append(paste(rep("SNP",num_SNPs),1:num_SNPs),"Class")




##############################################
# Specify pmin values
##############################################

pmin_list=getPminList()




#########################################################
# Start of Algorithm
#########################################################
timing_summary_table_CPU=matrix(ncol=2)
colnames(timing_summary_table_CPU)=c("Code Section","Time in Seconds")
timing_summary_table_system=matrix(ncol=2)
colnames(timing_summary_table_system)=c("Code Section","Time in Seconds")
timing_summary_table_user=matrix(ncol=2)
colnames(timing_summary_table_user)=c("Code Section","Time in Seconds")

diff_areas_alg_start_time <- proc.time()          #TIME

####################################################################################
# compute the p-values of every allele for all the classes
####################################################################################
computeSnpPvalue<-function(snp,classes)
{
  snp_col=current_fold[,snp]
  p_value=chisq.test(table(snp_col, classes))$p.value
  p_value=abs(log10(p_value))
  return(p_value)
}

generating_pvalues_tables_start_time <- proc.time()          #TIME

p_values_list=mapply(computeSnpPvalue,1:num_SNPs,MoreArgs=list(class_column))

generating_pvalues_tables_end_time <- proc.time()          #TIME
generating_pvalues_tables_duration <- generating_pvalues_tables_end_time-generating_pvalues_tables_start_time          #TIME

generating_pvalues_tables_duration_CPU=generating_pvalues_tables_duration['user.self']+generating_pvalues_tables_duration['sys.self']
generating_pvalues_tables_duration_CPU=round(generating_pvalues_tables_duration_CPU,digits=3)
generating_pvalues_tables_duration_CPU_row=c("Duration to Generate Pvalues Table",generating_pvalues_tables_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,generating_pvalues_tables_duration_CPU_row)

#TODO
# generating_pvalues_tables_duration_user=round(generating_pvalues_tables_duration['user.self'],digits=3)
# generating_pvalues_tables_duration_system=round(generating_pvalues_tables_duration['sys.self'],digits=3)
# generating_pvalues_tables_duration_elapsed=round(generating_pvalues_tables_duration['elapsed'],digits=3)



#########################################
# Identify Differentiable Locations
#########################################

identify_diff_locations_start_time <- proc.time()          #TIME

scaled_threshold_pvalue=abs(log10(0.05/num_SNPs))
all_dif_locations=which(p_values_list>=scaled_threshold_pvalue)

identify_diff_locations_end_time <- proc.time()          #TIME
identify_diff_locations_duration <- identify_diff_locations_end_time-identify_diff_locations_start_time          #TIME

identify_diff_locations_duration_CPU=identify_diff_locations_duration['user.self']+identify_diff_locations_duration['sys.self']
identify_diff_locations_duration_CPU=round(identify_diff_locations_duration_CPU,digits=3)
identify_diff_locations_duration_CPU_row=c("Duration to Identifying Differentiable Locations",identify_diff_locations_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,identify_diff_locations_duration_CPU_row)

#TODO
# identify_diff_locations_duration_user=round(identify_diff_locations_duration['user.self'],digits=3)
# identify_diff_locations_duration_system=round(identify_diff_locations_duration['sys.self'],digits=3)
# identify_diff_locations_duration_elapsed=round(identify_diff_locations_duration['elapsed'],digits=3)

########################################
# Creat list of differentiable areas
########################################

allowed_missing_SNPs=5 #(like window values...)

identify_diff_areas_start_time <- proc.time()          #TIME

differences=diff(all_dif_locations)
differentiable_areas=list()
i=1
while(i < length(all_dif_locations))
{
  start=i
  while((differences[i]<=allowed_missing_SNPs)&&(i < length(all_dif_locations)))
  {
    i=i+1
  }
  end=i
  motif=all_dif_locations[start:end]
  differentiable_areas[[length(differentiable_areas)+1]]=motif
  i=i+1
}

identify_diff_areas_end_time <- proc.time()          #TIME

identify_diff_areas_duration <- identify_diff_areas_end_time-identify_diff_areas_start_time          #TIME
identify_diff_areas_duration_CPU=identify_diff_areas_duration['user.self']+identify_diff_areas_duration['sys.self']
identify_diff_areas_duration_CPU=round(identify_diff_areas_duration_CPU,digits=3)
identify_diff_areas_duration_CPU_row=c("Duration to Identifying Differentiable Areas",identify_diff_areas_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,identify_diff_areas_duration_CPU_row)

#TODO
# identify_diff_areas_duration_user=round(identify_diff_areas_duration['user.self'],digits=3)
# identify_diff_areas_duration_system=round(identify_diff_areas_duration['sys.self'],digits=3)
# identify_diff_areas_duration_elapsed=round(identify_diff_areas_duration['elapsed'],digits=3)

######################
#Computing duration
######################

diff_areas_alg_end_time <- proc.time()

diff_areas_alg_time_table=diff_areas_alg_end_time-diff_areas_alg_start_time
diff_areas_alg_duration_CPU=diff_areas_alg_time_table['user.self']+diff_areas_alg_time_table['sys.self']
diff_areas_alg_duration_CPU=round(diff_areas_alg_duration_CPU,digits=3)
diff_areas_alg_duration_row=c("Total Duration of Differentiable Areas Algorithm",diff_areas_alg_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,diff_areas_alg_duration_row)

#TODO
# diff_areas_alg_duration_user=round(diff_areas_alg_time_table['user.self'],digits=3)
# diff_areas_alg_duration_system=round(diff_areas_alg_time_table['sys.self'],digits=3)
# diff_areas_alg_duration_elapsed=round(diff_areas_alg_time_table['elapsed'],digits=3)










##############################################################################################
# Get data for every area, generate pstree and identify motifs and their corresponding classes
##############################################################################################

#For building the trees
source(paste(home_dir,"/AlgorithmSrc/pstreeModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/cprobModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/maxRelN.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/mergeCprobModified.R",sep=""))

#For plotting the trees
source(paste(home_dir,"/AlgorithmSrc/plot_class_color_pstf.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plot_class_color_pstr.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plotNodeModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plotNodeProbModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plotPSTlegend.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plotPSTsetlayout.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/plotTreeModified.R",sep=""))
# source(paste(home_dir,"/AlgorithmSrc/as.pstree.updated.R",sep=""))

#For pruning the tree/get longest motif
source(paste(home_dir,"/AlgorithmSrc/getNodestoKeep.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/nodeKeepModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/isLeafSegmented.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/pruneUpdated.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/getPstrNode.R",sep=""))

list_all_motifs_form_one_diff_area_duration_CPU=c()
list_of_all_diff_motifs_alg_duration_CPU=c()
list_all_seq_def_duration_CPU=c()
list_all_pstree_modified_duration_CPU=c()
list_all_get_node_to_keep_duration_CPU=c()
list_all_prune_updated_duration_CPU=c()
list_all_motif_freq_positions_duration_CPU=c()
list_all_removing_overlap_motifs_duration_CPU=c()
list_of_all_predicting_all_diff_motifs_duration_CPU=c()

window=1

for(p in 1:length(pmin_list))
{
  
  #Timing
  diff_motifs_alg_start_time <- proc.time()
  predicting_all_diff_motifs_start_time <- proc.time()
  
  pmin=pmin_list[p]
  pmin_text=as.character(pmin)
  if(nchar(pmin_text)==3)
    pmin_text=paste(pmin_text,"0",sep="")

  all_predicted_motifs_table=matrix(ncol=5)
  motif_details_table_no_overlap=matrix(ncol=5)
  final_predicted_motifs_table=matrix(ncol=5)
  for(i in 1:length(differentiable_areas))  #For every differentiable area (TODO: Parallelize)
  {
    extract_motifs_form_one_diff_area_start_time <- proc.time() #Timing
    
    current_area=differentiable_areas[[i]]
    start=current_area[1]
    end=current_area[length(current_area)]
    if(start-window>=1)
      start=start-window
    if(end+window<=num_SNPs)
      end=end+window
    current_area=start:end
    current_area_label=paste(start,":",end,sep="")
    current_area_data=current_fold_no_class_column[current_area]
    end_column=rep("$",num_sequences)
    current_area_data=cbind(current_area_data,end_column)
    
    seq_def_start_time <- proc.time() #Timing
    current_area_data=seqdef(current_area_data)
    seq_def_end_time <- proc.time() #Timing
    
    seq_def_duration=seq_def_end_time-seq_def_start_time
    seq_def_duration_CPU=seq_def_duration['user.self']+seq_def_duration['sys.self']
    seq_def_duration_CPU=round(seq_def_duration_CPU,digits=3)
    list_all_seq_def_duration_CPU=c(list_all_seq_def_duration_CPU,seq_def_duration_CPU)
    
    pstree_modified_start_time <- proc.time() #Timing
    current_area_pst=pstreeModified(current_area_data,group=class_column,pmin=pmin)  
    pstree_modified_end_time <- proc.time() #Timing
    
    pstree_modified_duration=pstree_modified_end_time-pstree_modified_start_time
    pstree_modified_duration_CPU=pstree_modified_duration['user.self']+pstree_modified_duration['sys.self']
    pstree_modified_duration_CPU=round(pstree_modified_duration_CPU,digits=3)
    list_all_pstree_modified_duration_CPU=c(list_all_pstree_modified_duration_CPU,pstree_modified_duration_CPU)
    
    all_nodes=nodenames(current_area_pst)
    
    get_node_to_keep_start_time <- proc.time() #Timing
    motifs=getNodestoKeep(all_nodes,current_area_pst)
    get_node_to_keep_end_time <- proc.time() #Timing
    
    get_node_to_keep_duration=get_node_to_keep_end_time-get_node_to_keep_start_time
    get_node_to_keep_duration_CPU=get_node_to_keep_duration['user.self']+get_node_to_keep_duration['sys.self']
    get_node_to_keep_duration_CPU=round(get_node_to_keep_duration_CPU,digits=3)
    list_all_get_node_to_keep_duration_CPU=c(list_all_get_node_to_keep_duration_CPU,get_node_to_keep_duration_CPU)

    if(OutputPlots)
    {
      prune_updated_start_time <- proc.time() #Timing
      pruned_pst=pruneUpdated(current_area_pst,keep=motifs)
      prune_updated_end_time <- proc.time() #Timing
      
      prune_updated_duration=prune_updated_end_time-prune_updated_start_time
      prune_updated_duration_CPU=prune_updated_duration['user.self']+prune_updated_duration['sys.self']
      prune_updated_duration_CPU=round(prune_updated_duration_CPU,digits=3)
      list_all_prune_updated_duration_CPU=c(list_all_prune_updated_duration_CPU,prune_updated_duration_CPU)
      
    }
      
    

  
    motif_freq_positions_start_time <- proc.time() #Timing
    for(j in 1:length(motifs))
    {
      motif=motifs[j]
      motif_node=getPstrNode(motif,current_area_pst)
      motif_node_details=motif_node@index
      class=motif_node_details[,'group']
      if(length(class)==1)
      {
        relative_frequence_per_position=motif_node_details[,'rel_n']
        positions=motif_node_details[,4:ncol(motif_node_details)]
        positions=as.vector(na.omit(positions))
        positions=paste(positions,collapse=",")
        motif_details_row=c(class,motif,positions,current_area_label,round(relative_frequence_per_position,digits=3))
        all_predicted_motifs_table=rbind(all_predicted_motifs_table,motif_details_row)
      }
    }
    
    motif_freq_positions_end_time <- proc.time() #Timing
    
    motif_freq_positions_duration=motif_freq_positions_end_time-motif_freq_positions_start_time
    motif_freq_positions_duration_CPU=motif_freq_positions_duration['user.self']+motif_freq_positions_duration['sys.self']
    motif_freq_positions_duration_CPU=round(motif_freq_positions_duration_CPU,digits=3)
    list_all_motif_freq_positions_duration_CPU=c(list_all_motif_freq_positions_duration_CPU,motif_freq_positions_duration_CPU)
    
    
    extract_motifs_form_one_diff_area_end_time <- proc.time()
    extract_motifs_form_one_diff_area_duration=extract_motifs_form_one_diff_area_end_time-extract_motifs_form_one_diff_area_start_time
    
    extract_motifs_form_one_diff_area_duration_CPU=extract_motifs_form_one_diff_area_duration['user.self']+extract_motifs_form_one_diff_area_duration['sys.self']
    extract_motifs_form_one_diff_area_duration_CPU=round(extract_motifs_form_one_diff_area_duration_CPU,digits=3)
    
    list_all_motifs_form_one_diff_area_duration_CPU=c(list_all_motifs_form_one_diff_area_duration_CPU,extract_motifs_form_one_diff_area_duration_CPU)
    
    #TODO
#     extract_motifs_form_one_diff_area_duration_user=round(extract_motifs_form_one_diff_area_duration['user.self'],digits=3)
#     extract_motifs_form_one_diff_area_duration_system=round(extract_motifs_form_one_diff_area_duration['sys.self'],digits=3)
#     extract_motifs_form_one_diff_area_duration_elapsed=round(extract_motifs_form_one_diff_area_duration['elapsed'],digits=3)
    
  }#END OF for(i in 1:length(differentiable_areas))  #For every differentiable area (TODO: Parallelize)

  predicting_all_diff_motifs_end_time <- proc.time()
  predicting_all_diff_motifs_duration_table=predicting_all_diff_motifs_end_time-predicting_all_diff_motifs_start_time
  predicting_all_diff_motifs_duration_CPU=predicting_all_diff_motifs_duration_table['user.self']+predicting_all_diff_motifs_duration_table['sys.self']
  predicting_all_diff_motifs_duration_CPU=round(predicting_all_diff_motifs_duration_CPU,digits=3)
  list_of_all_predicting_all_diff_motifs_duration_CPU=c(list_of_all_predicting_all_diff_motifs_duration_CPU,predicting_all_diff_motifs_duration_CPU)
  
  #TODO
#   predicting_all_diff_motifs_duration_user=round(predicting_all_diff_motifs_duration_table['user.self'],digits=3)
#   predicting_all_diff_motifs_duration_system=round(predicting_all_diff_motifs_duration_table['sys.self'],digits=3)
#   predicting_all_diff_motifs_duration_elapsed=round(predicting_all_diff_motifs_duration_table['elapsed'],digits=3)    
    
    colnames(all_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
    rownames(all_predicted_motifs_table)=NULL

if(nrow(all_predicted_motifs_table)!=1)   #if motifs were found continue processing, else, output empty matrix
{
  if(nrow(all_predicted_motifs_table)==2)   #if only 1 motif was found
  {
    all_predicted_motifs_table=matrix(all_predicted_motifs_table[2,],nrow=1,ncol=5) #create a 1-row matrix

  }else if(nrow(all_predicted_motifs_table)>2)  #if more than 1 motif were found. remove NA row and order per class
  {
    all_predicted_motifs_table=all_predicted_motifs_table[-1,]
    all_predicted_motifs_table=all_predicted_motifs_table[order(all_predicted_motifs_table[,1]),]

  }

#********************************************************************************************************

removing_overlap_motifs_start_time <- proc.time() #Timing

class_column=all_predicted_motifs_table[,1]
motifs_column=all_predicted_motifs_table[,2]
positions_column=all_predicted_motifs_table[,3]
diff_areas_column=all_predicted_motifs_table[,4]
frequency_column=all_predicted_motifs_table[,5]
distinct_diff_areas=levels(factor(diff_areas_column))

for(i in 1:length(distinct_diff_areas))
{
  current_area=distinct_diff_areas[i]
  current_area_index=which(diff_areas_column==current_area)
  current_area_classes=class_column[current_area_index]
  distinct_current_area_classes=levels(factor(current_area_classes))
  for(n in 1:length(distinct_current_area_classes))
  {
    current_class=distinct_current_area_classes[n]
    same_area_same_class_motifs_index=which((diff_areas_column==current_area)&(class_column==current_class))
    if(length(same_area_same_class_motifs_index)==1)
    {
      current_row=all_predicted_motifs_table[same_area_same_class_motifs_index,]
      final_predicted_motifs_table=rbind(final_predicted_motifs_table,current_row)
    }
    else
    {
      motif_details_table_no_overlap=matrix(ncol=5)
      same_area_same_class_motifs=motifs_column[same_area_same_class_motifs_index]
      same_area_same_class_freq=frequency_column[same_area_same_class_motifs_index]
      same_area_same_class_positions=positions_column[same_area_same_class_motifs_index]
      
      max_length=0
      max_freq=0
      max_motif=""
      max_motif_start_position=0
      max_motif_end_location=0
      max_motif_details_row=0
      for(l in 1:length(same_area_same_class_motifs))
      {
        current_motif=same_area_same_class_motifs[l]
        current_motif_freq=as.numeric(same_area_same_class_freq[l])
        current_motif_length=nchar(current_motif)
        current_motif_start_position=as.numeric(same_area_same_class_positions[l])
        current_motif_end_position=current_motif_start_position+ceiling(current_motif_length/2)-1
        
        if(current_motif_length>max_length)
        {
          max_length=current_motif_length
          max_freq=current_motif_freq
          max_motif=current_motif
          max_motif_start_position=current_motif_start_position
          max_motif_end_position=current_motif_end_position
          max_motif_details_row_index=which((diff_areas_column==current_area)&(class_column==current_class)&(motifs_column==current_motif))
          max_motif_details_row=all_predicted_motifs_table[max_motif_details_row_index,]        
        }
        else if(current_motif_length==max_length)
        {
          #if length of current and max are equal, Check if there positions overlap:
          overlap=TRUE
          intersection=c()
          if((!is.na(max_motif_start_position))&&(!is.na(current_motif_start_position)))#to check if their are more than 1 start locations, then don't check for overlap
          {
            current_range=current_motif_start_position:current_motif_end_position
            max_range=max_motif_start_position:max_motif_end_position
            intersection=intersect(current_range,max_range)
            if(length(intersection)==0)
            {
              #keep both motifs
              overlap=FALSE
              #then, we keep track of the motif that will be removed
            }
          }
          if(current_motif_freq>max_freq) #remove max_motif
          {
            if(overlap==FALSE)  #if they do not overlap, store the motif that will be removed: max_motif in this case
            {
              no_overlap_row=max_motif_details_row
              motif_details_table_no_overlap=rbind(motif_details_table_no_overlap,no_overlap_row)
            }
            max_length=current_motif_length
            max_freq=current_motif_freq
            max_motif=current_motif
            max_motif_details_row_index=which((diff_areas_column==current_area)&(class_column==current_class)&(motifs_column==current_motif))
            max_motif_details_row=all_predicted_motifs_table[max_motif_details_row_index,]
          }
          else #remove current motif
          {
            if(overlap==FALSE)  #if they do not overlap, store the motif that will be removed: current_motif in this case
            {
              current_motif_details_row_index=which((diff_areas_column==current_area)&(class_column==current_class)&(motifs_column==current_motif))
              current_motif_details_row=all_predicted_motifs_table[current_motif_details_row_index,]
              no_overlap_row=current_motif_details_row
              motif_details_table_no_overlap=rbind(motif_details_table_no_overlap,no_overlap_row)
            }
          }
        }#END OF else if(current_motif_length==max_length)
      }# END OF for(l in 1:length(same_area_same_class_motifs)) LOOP
      
      final_predicted_motifs_table=rbind(final_predicted_motifs_table,max_motif_details_row)
      
      
      
      
      
      #Check if any of the removed motifs does not overlap with same class same area longer motifs
      for(o in 1:nrow(motif_details_table_no_overlap))
      {
        current_row=motif_details_table_no_overlap[o,]
        current_motif=current_row[2]
        if(!is.na(current_motif))
        {
          length_current_motif=nchar(current_motif)
          if(length_current_motif>=max_length)
          {
            final_predicted_motifs_table=rbind(final_predicted_motifs_table,current_row)
          }
        }
      }
      
      
    }#  END OF ELSE OF if(length(same_area_same_class_motifs_index)==1)
  }# END OF for(n in 1:length(distinct_current_area_classes)) LOOP
}# END of for(i in 1:length(distinct_diff_areas)) LOOP

removing_overlap_motifs_end_time <- proc.time() #Timing

removing_overlap_motifs_duration=removing_overlap_motifs_end_time-removing_overlap_motifs_start_time
removing_overlap_motifs_duration_CPU=removing_overlap_motifs_duration['user.self']+removing_overlap_motifs_duration['sys.self']
removing_overlap_motifs_duration_CPU=round(removing_overlap_motifs_duration_CPU,digits=3)
list_all_removing_overlap_motifs_duration_CPU=c(list_all_removing_overlap_motifs_duration_CPU,removing_overlap_motifs_duration_CPU)

colnames(final_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
rownames(final_predicted_motifs_table)=NULL

if(nrow(final_predicted_motifs_table)!=1)   #if motifs were found continue processing, else, output empty matrix
{
  if(nrow(final_predicted_motifs_table)==2)   #if only 1 motif was found
  {
    final_predicted_motifs_table=matrix(final_predicted_motifs_table[2,],nrow=1,ncol=5) #create a 1-row matrix
  }else if(nrow(final_predicted_motifs_table)>2)  #if more than 1 motif were found. remove NA row and order per class
  {
    final_predicted_motifs_table=final_predicted_motifs_table[-1,]
    final_predicted_motifs_table=final_predicted_motifs_table[order(final_predicted_motifs_table[,1]),]
  }
}


colnames(motif_details_table_no_overlap)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
rownames(motif_details_table_no_overlap)=NULL
motif_details_table_no_overlap=motif_details_table_no_overlap[order(motif_details_table_no_overlap[,1]),]  





#********************************************************************************************************


diff_motifs_alg_end_time <- proc.time()
diff_motifs_alg_duration_table=diff_motifs_alg_end_time-diff_motifs_alg_start_time
diff_motifs_alg_duration_CPU=diff_motifs_alg_duration_table['user.self']+diff_motifs_alg_duration_table['sys.self']
diff_motifs_alg_duration_CPU=round(diff_motifs_alg_duration_CPU,digits=3)
list_of_all_diff_motifs_alg_duration_CPU=c(list_of_all_diff_motifs_alg_duration_CPU,diff_motifs_alg_duration_CPU)

#TODO
#   diff_motifs_alg_duration_user=round(diff_motifs_alg_duration_table['user.self'],digits=3)
#   diff_motifs_alg_duration_system=round(diff_motifs_alg_duration_table['sys.self'],digits=3)
#   diff_motifs_alg_duration_elapsed=round(diff_motifs_alg_duration_table['elapsed'],digits=3)
}else{
  final_predicted_motifs_table=all_predicted_motifs_table
}

}#END OF PMIN_LIST LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


extract_motifs_form_one_diff_area_duration_row=c("Duration of Extracting Motifs from 1 Differentiable Area",round(mean(list_all_motifs_form_one_diff_area_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,extract_motifs_form_one_diff_area_duration_row)



seq_def_duration_row=c("Duration of seqdef",round(mean(list_all_seq_def_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,seq_def_duration_row)


pstree_modified_duration_row=c("Duration of pstreeModified",round(mean(list_all_pstree_modified_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,pstree_modified_duration_row)

get_node_to_keep_duration_row=c("Duration of getNodeToKeep",round(mean(list_all_get_node_to_keep_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,get_node_to_keep_duration_row)

if(OutputPlots)
{
  prune_updated_duration_row=c("Duration of pruneUpdated",round(mean(list_all_prune_updated_duration_CPU),digits=3))
  timing_summary_table_CPU=rbind(timing_summary_table_CPU,prune_updated_duration_row)
}


motif_freq_positions_duration_row=c("Duration of Identifying Motif Locations and Frequencies",round(mean(list_all_motif_freq_positions_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,motif_freq_positions_duration_row)



predicting_all_diff_motifs_duration_row=c("Duration of Predicting all Diff Motifs",round(mean(list_of_all_predicting_all_diff_motifs_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,predicting_all_diff_motifs_duration_row)

removing_overlap_motifs_duration_row=c("Duration of Removing Overlapping Motifs",round(mean(list_all_removing_overlap_motifs_duration_CPU),digits=3))
timing_summary_table_CPU=rbind(timing_summary_table_CPU,removing_overlap_motifs_duration_row)


average_diff_motifs_alg_duration_CPU=round(mean(list_of_all_diff_motifs_alg_duration_CPU),digits=3)
average_diff_motifs_alg_duration_row=c("Duration of Whole Differentiable Motifs Algorithm",average_diff_motifs_alg_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,average_diff_motifs_alg_duration_row)


total_alg_duration_CPU=round(sum(diff_areas_alg_duration_CPU,average_diff_motifs_alg_duration_CPU),digits=3)

total_alg_duration_row=c("Total Algorithm Duration",total_alg_duration_CPU)
timing_summary_table_CPU=rbind(timing_summary_table_CPU,total_alg_duration_row)

timing_summary_table_CPU=timing_summary_table_CPU[-1,]

# current_output_file_name_short=getOutputFileNameShort()
print("current_output_file_name_short")
print(current_output_file_name_short)

time_file_name=""
if(OutputRuns)
{
  time_file_name=paste(getRunTimeEvaluationPerFoldDir(),"/AlgorithmRuntime_",current_output_file_name_short,".csv",sep="")
  if(data_type=="complete")
  {
    time_file_name=paste(getRunTimeEvaluationDir(),"/AlgorithmRuntime_",current_output_file_name_short,".csv",sep="")
  }
}else{
  time_file_name=paste(getRunTimeEvaluationDir(),"/AlgorithmRuntime_",current_output_file_name_short,".csv",sep="")
}
write.csv(timing_summary_table_CPU,time_file_name,row.names=FALSE)

print("ALGORITHM CODE for timing ENDED!") 

time=proc.time() - ptm  
print(paste("TIME in sec: ",time['elapsed']))

