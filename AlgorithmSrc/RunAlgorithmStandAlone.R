###########################################
# RunAlgorithm
# Last updated: 23-02-15
###########################################

print("inside RunAlgorithmStandAlone.R")

library("ggplot2")
library("TraMineR") 
library("PST") 
library("plyr")

ptm <- proc.time()

source_code_dir=getwd()   #directory of the algorithm's source code
args <- commandArgs(trailingOnly = TRUE)
file_name=args[1]
pmin_list=args[2]
pmin_list=unlist(strsplit(pmin_list,","))
output_dir=args[3]
verbose=args[4]    #Gioia
data_set_name=basename(file_name)

diff_locations_dir=file.path(output_dir,"6-PredictedDifferentiableAreas")
dir.create(diff_locations_dir, showWarnings = FALSE)

trees_dir=file.path(output_dir,"10-Trees")
dir.create(trees_dir, showWarnings = FALSE)

predicted_motifs_dir=file.path(output_dir,"8-PredictedMotifs")
dir.create(predicted_motifs_dir, showWarnings = FALSE)

#############################################################
#Import data and then fix column names:
#############################################################

data <- read.table(file_name, quote="\"")
# data <- read.table(file_name, quote="\"",header=TRUE)


class_column=data[,ncol(data)]
data_no_class_column=data[,-length(data)]
colnames(data_no_class_column)=1:ncol(data_no_class_column)
num_sequences= nrow(data)             
num_SNPs= ncol(data)-1    


############################################
# Compute all the p-values
############################################

print("computing p-values")

computeSnpPvalue<-function(snp,classes)
{
#   print("snp")
#   print(snp)
  snp_col=data[,snp]
#   print("snp_col")
#   print(snp_col)
#   print("table(snp_col,classes)")
#   print(table(snp_col,classes))
  p_value=chisq.test(table(snp_col, classes))$p.value
#   print("p_value")
#   print(p_value)
  p_value=abs(log10(p_value))
#   print("abs(p_value)")
#   print(p_value)
  return(p_value)
}



p_values_list=mapply(computeSnpPvalue,1:num_SNPs,MoreArgs=list(class_column))
# print(warnings())
# print("p_values_list")
# print(p_values_list)
# print("max(p_values_list)")
# print(max(p_values_list))



#########################################
# Identify Differentiable Locations
#########################################
print("Extracting Differentiable locations")
scaled_threshold_pvalue=abs(log10(0.05/num_SNPs))
# print("scaled_threshold_pvalue")
# print(scaled_threshold_pvalue)
all_dif_locations=which(p_values_list>=scaled_threshold_pvalue)
# print("all_dif_locations")
# print(all_dif_locations)

if(length(all_dif_locations)==0)
{
  print("There were no differentiable locations found.")
  quit()
}

if(verbose)
{
  
  write(all_dif_locations,paste(diff_locations_dir,"/DifferentiableLocations_",data_set_name,".txt",sep=""))
}

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

write.table(differentiable_areas_matrix,paste(diff_locations_dir,"/DifferentiableAreas_",data_set_name,".csv",sep=""), sep=",",na=""
            ,col.names=FALSE,row.names=FALSE)


##############################################################################################
# Get data for every area, generate pstree and identify motifs and their corresponding classes
##############################################################################################

#For building the trees
source(paste(source_code_dir,"/pstreeModified.R",sep=""))
source(paste(source_code_dir,"/cprobModified.R",sep=""))
source(paste(source_code_dir,"/maxRelN.R",sep=""))
source(paste(source_code_dir,"/mergeCprobModified.R",sep=""))

#For plotting the trees
source(paste(source_code_dir,"/plot_class_color_pstf.R",sep=""))
source(paste(source_code_dir,"/plot_class_color_pstr.R",sep=""))
source(paste(source_code_dir,"/plotNodeModified.R",sep=""))
source(paste(source_code_dir,"/plotNodeProbModified.R",sep=""))
source(paste(source_code_dir,"/plotPSTlegend.R",sep=""))
source(paste(source_code_dir,"/plotPSTsetlayout.R",sep=""))
source(paste(source_code_dir,"/plotTreeModified.R",sep=""))

#For pruning the tree/get longest motif
source(paste(source_code_dir,"/getNodestoKeep.R",sep=""))
source(paste(source_code_dir,"/nodeKeepModified.R",sep=""))
source(paste(source_code_dir,"/isLeafSegmented.R",sep=""))
source(paste(source_code_dir,"/pruneUpdated.R",sep=""))
source(paste(source_code_dir,"/getPstrNode.R",sep=""))

print("Extracting Differentiable Motifs")


window=1
min_motif_length=1

for(p in 1:length(pmin_list))
{
  pmin=pmin_list[p]
  pmin_text=as.character(pmin)
  if(nchar(pmin_text)==3)
    pmin_text=paste(pmin_text,"0",sep="")
  
  trees_file=paste(trees_dir,"/Tree_",data_set_name,"_pmin:",pmin_text,".pdf",sep="")
  all_predicted_motifs_file=paste(predicted_motifs_dir,"/AllPredictedMotifs_",data_set_name,"_pmin:",pmin_text,".csv",sep="")
  final_predicted_motifs_file=paste(predicted_motifs_dir,"/PredictedMotifs_",data_set_name,"_pmin:",pmin_text,".csv",sep="")
  motif_details_table_no_overlap_file=paste(predicted_motifs_dir,"/NoOverlapMotifs_",data_set_name,"_pmin:",pmin_text,".csv",sep="")
  pdf(trees_file) 
  
  all_predicted_motifs_table=matrix(ncol=5)
  motif_details_table_no_overlap=matrix(ncol=5)
  final_predicted_motifs_table=matrix(ncol=5)
  for(i in 1:length(differentiable_areas))  #For every differentiable area (TODO: Parallelize)
  {
#     print(paste("Differentiable area",i))
    current_area=differentiable_areas[[i]]
#     print("current_area")
#     print(current_area)
    if(length(current_area)>=min_motif_length)
    {
      start=current_area[1]
      end=current_area[length(current_area)]
      if(start-window>=1)
        start=start-window
      if(end+window<=num_SNPs)
        end=end+window
      current_area=start:end
      current_area_label=paste(start,":",end,sep="")
      current_area_data=data_no_class_column[current_area]
      end_column=rep("$",num_sequences)
      current_area_data=cbind(current_area_data,end_column)
#       print("current_area_data")
#       print(head(current_area_data))
#       print("Creating the seqdata object")
      current_area_data=seqdef(current_area_data)
      #     print("current_area_data after seqdef")
      #     print(current_area_data)
#       print("Generating PST")
      current_area_pst=pstreeModified(current_area_data,group=class_column,pmin=pmin)    #TODO: try various pmins
      all_nodes=nodenames(current_area_pst)
#       print("current area PST nodes")
#       print(all_nodes)
      if(length(all_nodes)!=1)    #If no differentiable motif was found for this differentiable area... go to next area
      {
        motifs=getNodestoKeep(all_nodes,current_area_pst)
#         print("nodes/motifs to keep:")
#         print(motifs)
#         print("prune tree and keep only the above nodes")
#         pruned_pst=pruneUpdated(current_area_pst,keep=motifs)
#         
#         if(verbose)
#         {
#           print("inside if verbose before plot_class_color_pstf")
#           plot_class_color_pstf(pruned_pst,title=paste("PRUNED",current_area_label))
#         }
#         print("Plot complete tree")
#         plot_class_color_pstf(current_area_pst,title=paste("COMPLETE",current_area_label))
        
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
        }#end of for(j in 1:length(motifs))
      }#end of length(all_nodes!=1)
      else
      {
        print("There were no diff. motifs found for the current differentiable area")
      }
    }#end of if(length(current_area)>=min_motif_length)
    else
    {
      print("The current area is very small. so do not search for motif, skip to next area...")
    }
  }# end of for(i in 1:length(differentiable_areas))
  
  dev.off (); 
  colnames(all_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
  rownames(all_predicted_motifs_table)=NULL
  all_predicted_motifs_table=all_predicted_motifs_table[-1,]
  all_predicted_motifs_table=all_predicted_motifs_table[order(all_predicted_motifs_table[,1]),]
  if(verbose)
  {
    write.csv(all_predicted_motifs_table,all_predicted_motifs_file,row.names=FALSE)
  }
  
#   print("nrow(all_predicted_motifs_table)")
#   print(nrow(all_predicted_motifs_table))
  
  if(nrow(all_predicted_motifs_table)!=0)   # if we do have differentiable motifs...
  {
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
    
    colnames(final_predicted_motifs_table)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
    rownames(final_predicted_motifs_table)=NULL
    final_predicted_motifs_table=final_predicted_motifs_table[-1,]
    final_predicted_motifs_table=final_predicted_motifs_table[order(final_predicted_motifs_table[,1]),]  
    write.csv(final_predicted_motifs_table,final_predicted_motifs_file,row.names=FALSE)
    
    
    if(verbose)
    {
      colnames(motif_details_table_no_overlap)=c("Class","Motif","Positions","Differentiable Area","Relative Frequency")
      rownames(motif_details_table_no_overlap)=NULL
      motif_details_table_no_overlap=motif_details_table_no_overlap[order(motif_details_table_no_overlap[,1]),]
      write.csv(motif_details_table_no_overlap,motif_details_table_no_overlap_file,row.names=FALSE)
    }
    
  }#end of if(nrow(all_predicted_motifs_table)!=0)  
  else
  {
    print("There were no differentiable motifs found. Only differentibale locations")
  }
  
}#END OF PMIN_LIST LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print("ALGORITHM CODE ENDED!") 

time=proc.time() - ptm  
print(paste("TIME in sec: ",time['elapsed']))



