###########################################
# RunAlgorithm
# Last updated: 23-02-15
###########################################

print("inside RunAlgorithm.R")

home_dir=getHomeDir()
OutputText=getOutputText()
OutputPlots=getOutputPlots()

ptm <- proc.time()


pmin_list=getPminListText()
pmin_list=unlist(strsplit(pmin_list,","))

output_dir=getAllOutputDir()

#############################################################
#Import data and then fix column names:
#############################################################
# data <- read.table(file_name, quote="\"")


# if(getDataType()=="UnseenTesting")
# {
#   current_fold_name="UnseenTesting"
# }else{
  current_fold_name=getCurrentFoldName()
# }

current_run_name=getCurrentRunName()

current_fold=getCurrentFold()

num_seq=nrow(current_fold)
print("num_seq")
print(num_seq)

class_column=current_fold[,ncol(current_fold)]
current_fold_no_class_column=current_fold[,-length(current_fold)]
colnames(current_fold_no_class_column)=1:ncol(current_fold_no_class_column)
num_sequences= nrow(current_fold)             
num_SNPs= ncol(current_fold)-1    

class_labels_column=current_fold[,ncol(current_fold)]
class_details=t(matrix(table(class_labels_column)))
distinct_classes=levels(as.factor(class_labels_column))
colnames(class_details)=distinct_classes

############################################
# Compute all the p-values
############################################

print("computing p-values")

computeSnpPvalue<-function(snp,classes)
{
  snp_col=current_fold[,snp]
  p_value=chisq.test(table(snp_col, classes))$p.value
  p_value=abs(log10(p_value))
  return(p_value)
}

print("current_output_file_name_short")
print(current_output_file_name_short)
# print("output_file_name_long")
# print(output_file_name_long)

p_values_list=mapply(computeSnpPvalue,1:num_SNPs,MoreArgs=list(class_column))
p_values_row=append(current_output_file_name_short,p_values_list)
p_values_matrix=matrix(p_values_row,nrow=1,byrow=TRUE)

file_pvalues_all_alleles=paste(getPvaluesDataPerAlleleDir(),"/PvaluesAllAlleles.txt",sep="")
write.table(p_values_matrix,file=file_pvalues_all_alleles,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)



#########################################
# Identify Differentiable Locations
#########################################
print("Extracting Differentiable locations")
scaled_threshold_pvalue=abs(log10(0.05/num_SNPs))
all_dif_locations=which(p_values_list>=scaled_threshold_pvalue)

if(OutputText)
{
  write(all_dif_locations,paste(getPredictedDifferentiableAreasDir(),"/DifferentiableLocations_",current_output_file_name_short,".txt",sep=""))
}

print("length(all_dif_locations)")
print(length(all_dif_locations))
########################################
# Creat list of differentiable areas
########################################
print("Extracting Differentiable areas")

allowed_missing_SNPs=5 #(like window values...)

differences=diff(all_dif_locations)
differentiable_areas=list()
differentiable_areas_matrix=matrix()
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
  motif_matrix=matrix(motif,nrow=1,ncol=length(motif))
  differentiable_areas_matrix=rbind.fill.matrix(differentiable_areas_matrix,motif_matrix)
  differentiable_areas[[length(differentiable_areas)+1]]=motif

  i=i+1
}

# lapply(differentiable_areas, write, paste(getPredictedDifferentiableAreasDir(),"/DifferentiableAreas_",current_output_file_name_short,".txt",sep=""), append=TRUE, ncolumns=20)
# lapply(differentiable_areas, write.table, paste(getPredictedDifferentiableAreasDir(),"/DifferentiableAreas_",current_output_file_name_short,".csv",sep=""), append=TRUE,sep=",")
write.table(differentiable_areas_matrix,paste(getPredictedDifferentiableAreasDir(),"/DifferentiableAreas_",current_output_file_name_short,".csv",sep=""), sep=",",na=""
            ,col.names=FALSE,row.names=FALSE)


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

#For pruning the tree/get longest motif
source(paste(home_dir,"/AlgorithmSrc/getNodestoKeep.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/nodeKeepModified.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/isLeafSegmented.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/pruneUpdated.R",sep=""))
source(paste(home_dir,"/AlgorithmSrc/getPstrNode.R",sep=""))
 
window=1

# getWindow<<-function()
# {
#   return(window)
# }
print("Running Algorithm")

for(p in 1:length(pmin_list))
{
  pmin=pmin_list[p]
  pmin_text=as.character(pmin)
  if(nchar(pmin_text)==3)
    pmin_text=paste(pmin_text,"0",sep="")

  trees_file=paste(getTreesDir(),"/Tree_",current_output_file_name_short,"_pmin_",pmin_text,".pdf",sep="")
  all_predicted_motifs_file=paste(getPredictedMotifsDir(),"/AllPredictedMotifs_",current_output_file_name_short,"_pmin_",pmin_text,".csv",sep="")
  final_predicted_motifs_file=paste(getPredictedMotifsDir(),"/PredictedMotifs_",current_output_file_name_short,"_pmin_",pmin_text,".csv",sep="")
  filtered_predicted_motifs_file=paste(getPredictedMotifsDir(),"/FilteredPredictedMotifs_",current_output_file_name_short,"_pmin_",pmin_text,".csv",sep="")
  filtered_predicted_motifs_all_folds_file=paste(getMotifsOutputDir(),"/",getDataType(),"/FilteredPredictedMotifsAllFolds_",current_output_file_name_short,"_pmin_",pmin_text,".csv",sep="")
  motif_details_table_no_overlap_file=paste(getPredictedMotifsDir(),"/NoOverlapMotifs_",current_output_file_name_short,"_pmin_",pmin_text,".csv",sep="")
  if(OutputRuns)
  {
    pdf(trees_file) 
  }

  all_predicted_motifs_table=matrix(ncol=5)
  motif_details_table_no_overlap=matrix(ncol=5)
  final_predicted_motifs_table=matrix(ncol=5)
  print("differentiable_areas")
  print(differentiable_areas)
  print("length(differentiable_areas)")
  print(length(differentiable_areas))
  if(length(differentiable_areas)>0)
  {
    for(i in 1:length(differentiable_areas))  #For every differentiable area (TODO: Parallelize)
    {
      current_area=differentiable_areas[[i]]
      start=current_area[1]
      end=current_area[length(current_area)]
      if(start-window>=1)
        start=start-window
      if(end+window<=num_SNPs)
        end=end+window
      current_area=start:end
      print("current_area")
      print(current_area)
      current_area_label=paste(start,":",end,sep="")
      current_area_data=current_fold_no_class_column[current_area]
      end_column=rep("$",num_sequences)
      current_area_data=cbind(current_area_data,end_column)
      current_area_data=seqdef(current_area_data)
      current_area_pst=pstreeModified(current_area_data,group=class_column,pmin=pmin)    #TODO: try various pmins
      
      all_nodes=nodenames(current_area_pst)
      print("all_nodes")
      print(all_nodes)
      motifs=getNodestoKeep(all_nodes,current_area_pst)
      print("motifs")
      print(motifs)
      
      if(OutputPlots)
      {
        pruned_pst=pruneUpdated(current_area_pst,keep=motifs)
        print("pruned_pst")
        print(pruned_pst)
        plot_class_color_pstf(pruned_pst,title=paste("PRUNED",current_area_label))
      }
      if(OutputRuns)
      {
        plot_class_color_pstf(current_area_pst,title=paste("COMPLETE",current_area_label))
      }
      
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
        else
        {
          print("the current motif is found in more that 1 class, so do not include in the table")
        }
      }
    }
    if(OutputRuns)
    {
      dev.off ()
    }
  }else{
    print("--------------------------------no diff areas!!")
    }#end of if(length(differentiable_areas)>0)
  

  print("all_predicted_motifs_table")
  print(all_predicted_motifs_table)
  
  colnames(all_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
  rownames(all_predicted_motifs_table)=NULL
  
  if(nrow(all_predicted_motifs_table)!=1)   #if motifs were found continue processing, else, output empty matrix
  {
    if(nrow(all_predicted_motifs_table)==2)   #if only 1 motif was found
    {
      all_predicted_motifs_table=matrix(all_predicted_motifs_table[2,],nrow=1,ncol=5) #create a 1-row matrix
      print("all_predicted_motifs_table")
      print(all_predicted_motifs_table)
      colnames(all_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
      if(OutputText)
      {
        write.csv(all_predicted_motifs_table,all_predicted_motifs_file,row.names=FALSE)
      }
    }else if(nrow(all_predicted_motifs_table)>2)  #if more than 1 motif were found. remove NA row and order per class
    {
      all_predicted_motifs_table=all_predicted_motifs_table[-1,]
      all_predicted_motifs_table=all_predicted_motifs_table[order(all_predicted_motifs_table[,1]),]
      if(OutputText)
      {
        write.csv(all_predicted_motifs_table,all_predicted_motifs_file,row.names=FALSE)
      }
    }
    
    
    #Continue processing after fixing matrix
    class_column=all_predicted_motifs_table[,1]
    motifs_column=all_predicted_motifs_table[,2]
    positions_column=all_predicted_motifs_table[,3]
    diff_areas_column=all_predicted_motifs_table[,4]
    frequency_column=all_predicted_motifs_table[,5]
    distinct_diff_areas=levels(factor(diff_areas_column))
    
    for(i in 1:length(distinct_diff_areas)) #for every differentiable area found
    {
      current_area=distinct_diff_areas[i]
      current_area_index=which(diff_areas_column==current_area)
      current_area_classes=class_column[current_area_index]
      distinct_current_area_classes=levels(factor(current_area_classes))
      for(n in 1:length(distinct_current_area_classes))
      {
        current_class=distinct_current_area_classes[n]
        same_area_same_class_motifs_index=which((diff_areas_column==current_area)&(class_column==current_class))
        if(length(same_area_same_class_motifs_index)==1)# if none of the motifs for this area have the same class lable, add to the final output
        {
          current_row=all_predicted_motifs_table[same_area_same_class_motifs_index,]
          final_predicted_motifs_table=rbind(final_predicted_motifs_table,current_row)
        }
        else
        {# there are 2 motifs for the same are and same class
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
            
            if(current_motif_length>max_length) # if they differ in length, keep the longest motif
            {
              max_length=current_motif_length
              max_freq=current_motif_freq
              max_motif=current_motif
              max_motif_start_position=current_motif_start_position
              max_motif_end_position=current_motif_end_position
              max_motif_details_row_index=which((diff_areas_column==current_area)&(class_column==current_class)&(motifs_column==current_motif))
              max_motif_details_row=all_predicted_motifs_table[max_motif_details_row_index,]        
            }
            else if(current_motif_length==max_length) # if they have the same length.
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
    
    print("final_predicted_motifs_table")
    print(final_predicted_motifs_table)
    colnames(final_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
    rownames(final_predicted_motifs_table)=NULL
    
    if(nrow(final_predicted_motifs_table)!=1)   #if motifs were found continue processing, else, output empty matrix
    {
      if(nrow(final_predicted_motifs_table)==2)   #if only 1 motif was found
      {
        final_predicted_motifs_table=matrix(final_predicted_motifs_table[2,],nrow=1,ncol=5) #create a 1-row matrix
        print("final_predicted_motifs_table")
        print(final_predicted_motifs_table)
        colnames(final_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
      }else if(nrow(final_predicted_motifs_table)>2)  #if more than 1 motif were found. remove NA row and order per class
      {
        final_predicted_motifs_table=final_predicted_motifs_table[-1,]
        final_predicted_motifs_table=final_predicted_motifs_table[order(final_predicted_motifs_table[,1]),]
      }
    }
    
    write.csv(final_predicted_motifs_table,final_predicted_motifs_file,row.names=FALSE)

    if(OutputText)
    {
      print("final_predicted_motifs_table")
      print(final_predicted_motifs_table)
      colnames(motif_details_table_no_overlap)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
      rownames(motif_details_table_no_overlap)=NULL
      motif_details_table_no_overlap=motif_details_table_no_overlap[order(motif_details_table_no_overlap[,1]),]
      write.csv(motif_details_table_no_overlap,motif_details_table_no_overlap_file,row.names=FALSE)
    }
    
  } else {
  final_predicted_motifs_table=all_predicted_motifs_table
  if(OutputText)
  {
    write.csv(all_predicted_motifs_table,all_predicted_motifs_file,row.names=FALSE)
  }
  write.csv(final_predicted_motifs_table,final_predicted_motifs_file,row.names=FALSE)
}
  
}#END OF PMIN_LIST LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















###############################################
# Computing frequency of motifs per class 
###############################################

upper_bound=getFoldDifferenceUpperBound()
lower_bound=getFoldDifferenceLowerBound()
motif_length_lower_bound=getMinMotifLength()



start_col=final_predicted_motifs_table[,3]
# end_col=final_predicted_motifs_table[,5]
motif_col=final_predicted_motifs_table[,2]
class_col=final_predicted_motifs_table[,1]


fold_change_col_name=c()
for (d in 1:(length(distinct_classes)-1))
{
  current_class=distinct_classes[d]
  print("current_class")
  print(current_class)
  for (b in (d+1):length(distinct_classes))
  {
    second_class=distinct_classes[b]
    print("second_class")
    print(second_class)
    corresponding_classes=paste(current_class,"div",second_class,sep="")
    print("corresponding_classes")
    print(corresponding_classes)
    fold_change_col_name=c(fold_change_col_name,corresponding_classes)
  }
}

num_fold_change_combinations=length(fold_change_col_name)

freq_table=matrix(ncol=(length(distinct_classes)+9+num_fold_change_combinations),nrow=1)   



for(i in 1:length(motif_col))
{
  current_motif=motif_col[i]
  if(!is.na(current_motif))
  {
    print("current_motif###################################################################")
    print(current_motif)
    current_motif_length=length(unlist(strsplit(current_motif,"-")))
    print("current_motif_length")
    print(current_motif_length)
    majority_class=class_col[i]
    print("majority_class")
    print(majority_class)
    cur_loc_start=as.character(start_col[i])
    print("cur_loc_start")
    print(cur_loc_start)
    split_loc_start=unlist(strsplit(cur_loc_start,","))
#     cur_loc_end=as.numeric(cur_loc_start)+length(split_loc_start)-1
#     print('cur_loc_end')
#     print(cur_loc_end)
#     cur_loc_end=as.character(end_col[i])
#     split_loc_end=unlist(strsplit(cur_loc_end,"-"))
    current_start_loc=0
    current_end_loc=0
    
    for(x in 1:length(split_loc_start))
    {
      print("location separator------------------------")
      current_start_loc=split_loc_start[x]
      print("current_start_loc")
      print(current_start_loc)
      current_end_loc=as.numeric(current_start_loc)+current_motif_length-1
      print("current_end_loc")
      print(current_end_loc)
      current_range=current_start_loc:current_end_loc
      #     current_range_rs_id=snp_rs_id[current_range]
      current_range_data=current_fold_no_class_column[,current_range]
      
      if(nchar(current_motif)==1)
      {
        print("inside if(nchar(current_motif)==1)")
        print("class(current_range_data)")
        print(class(current_range_data))
        current_range_data=as.matrix(current_range_data)
        #       print("current_range_data")
        #       print(head(current_range_data))
      }
      
      print("dim(current_range_data)")
      print(dim(current_range_data))
      
      all_classes_freqs=c()
      print("distinct_classes")
      print(distinct_classes)
      for(l in 1:length(distinct_classes))
      {
        current_class=distinct_classes[l]
        print("current_class")
        print(current_class)
        print("class(current_class)")
        print(class(current_class))
        print("length(class_labels_column)")
        print(length(class_labels_column))
        print("table(class_labels_column)")
        print(table(class_labels_column))
        current_class_indexes=which(class_labels_column==current_class)
        print("length(current_class_indexes)")
        print(length(current_class_indexes))
        current_class_range_data=current_range_data[current_class_indexes,]
        print("dim(current_class_range_data)")
        print(dim(current_class_range_data))
        
        if(nchar(current_motif)==1)
        {
          print("inside if(nchar(current_motif)==1)")
          print("class(current_class_range_data)")
          print(class(current_class_range_data))
          current_class_range_data=as.matrix(current_class_range_data)
          #         print("current_class_range_data")
          #         print(head(current_class_range_data))
        }
        
        print("dim(current_class_range_data)")
        print(dim(current_class_range_data))
        
        current_class_num_seq=nrow(current_class_range_data)
        print("current_class_num_seq")
        print(current_class_num_seq)
        
        current_class_motifs=apply(current_class_range_data,1,paste,collapse="-")
        print("length(current_class_motifs)")
        print(length(current_class_motifs))
        
        current_class_motifs_freqs=table(current_class_motifs)
        print("current_class_motifs_freqs")
        print(current_class_motifs_freqs)
        current_motif_freq_per_class=current_class_motifs_freqs[current_motif]
        print("current_motif_freq_per_class")
        print(current_motif_freq_per_class)
        current_motif_rel_freq_per_class=round(current_motif_freq_per_class/current_class_num_seq,2)
        print("current_motif_rel_freq_per_class")
        print(current_motif_rel_freq_per_class)
        if(is.na(current_motif_rel_freq_per_class))
        {
          current_motif_rel_freq_per_class=0
        }
        all_classes_freqs=c(all_classes_freqs,current_motif_rel_freq_per_class)
      }#END of for(l in 1:length(distinct_classes))
      
      print("current_fold_name")
      print(current_fold_name)
      print("current_motif")
      print(current_motif)
      current_motif_length=ceiling(nchar(current_motif)/2)
      print("current_motif_length")
      print(current_motif_length)
      print("majority_class")
      print(majority_class)
      print("current_start_loc")
      print(current_start_loc)
      print("current_end_loc")
      print(current_end_loc)
      print("all_classes_freqs")
      print(all_classes_freqs)
      fold_change_vec=c()
      for (d in 1:(length(all_classes_freqs)-1))
      {
        current_class_freq=as.numeric(all_classes_freqs[d])
        print("current_class_freq")
        print(current_class_freq)
        #       current_class=distinct_classes[d]
        #       print("current_class")
        #       print(current_class)
        for (b in (d+1):length(all_classes_freqs))
        {
          second_class_freq=as.numeric(all_classes_freqs[b])
          print("second_class_freq")
          print(second_class_freq)
          #         second_class=distinct_classes[b]
          #         print("second_class")
          #         print(second_class)
          fold_change=current_class_freq/second_class_freq
          if((second_class_freq==0)&&(current_class_freq==0))
          {
            print("Numerator and Dinominator are 0")
            fold_change=-1
          }else if(second_class_freq==0){
            print("Donominator only = 0")
            fold_change=999999
            
          }
          print("fold_change")
          print(fold_change)
          fold_change_vec=c(fold_change_vec,fold_change)
          #         corresponding_classes=paste(current_class,"div",second_class,sep="")
          #         print("corresponding_classes")
          #         print(corresponding_classes)
          #         fold_change_col_name=c(fold_change_col_name,corresponding_classes)
        }#end of for (b in (d+1):length(all_classes_freqs))
      }# end of for (d in 1:(length(all_classes_freqs)-1))
      
      print("fold_change_vec")
      print(fold_change_vec)
      print("fold_change_col_name")
      print(fold_change_col_name)
      max_fold_change=max(fold_change_vec)
      print("max_fold_change")
      print(max_fold_change)
      min_fold_change=min(fold_change_vec)
      print("min_fold_change")
      print(min_fold_change)
      freq_table_row=c(current_run_name,current_fold_name,majority_class,current_motif,current_motif_length,
                       current_start_loc,current_end_loc,all_classes_freqs,
                       fold_change_vec,max_fold_change,min_fold_change)
#       if(((max_fold_change>upper_bound)||(min_fold_change<lower_bound))&&(as.numeric(current_motif_length)>=as.numeric(motif_length_lower_bound)))
#         {
        freq_table=rbind(freq_table,freq_table_row)
        print("row added")
#         }
      }#END of for(x in 1:length(split_loc_start))
    }else{
    freq_table_row=c(current_run_name,current_fold_name,rep("__",5),rep("__",length(distinct_classes)),rep("__",length(fold_change_col_name)),"__","__")
    freq_table=rbind(freq_table,freq_table_row)
      }# End of if(!is.na(current_motif))
}








# if(nrow(freq_table)==1)   # that is no motif was found
# {
#   freq_table_row=c(current_fold_name,"__","__","__","__","__","__","__","__","__","__")
#   freq_table=rbind(freq_table,freq_table_row)
# }

freq_table=freq_table[-1,]

if(class(freq_table)!='matrix')
{
  freq_table=t(as.matrix(freq_table))
  print("freq_table 1 row")
  print(freq_table)
}else{
  print("freq_table")
  print(head(freq_table))
  freq_table=freq_table[order(as.numeric(freq_table[,5])),]
  print("freq_table after order")
  print(freq_table)
}



colnames(freq_table)=c("Run","Fold","Majority_Class","Motif","Motif_Length","Start_Location",
                       "End_Location",paste("Freq_Class_",distinct_classes,sep=""),
                       fold_change_col_name,"Max_Fold_Change","Min_Fold_Change")

# write.csv(freq_table,filtered_predicted_motifs_file,row.names=F,col.names=T)

if(file.exists(filtered_predicted_motifs_all_folds_file)==FALSE) {
  write.table(freq_table,file=filtered_predicted_motifs_all_folds_file,sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE,append=TRUE)
} else if(file.exists(filtered_predicted_motifs_all_folds_file)==TRUE) {
  write.table(freq_table,file=filtered_predicted_motifs_all_folds_file,sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  
}












print("ALGORITHM CODE ENDED!") 

time=proc.time() - ptm  
print(paste("TIME in sec: ",time['elapsed']))




