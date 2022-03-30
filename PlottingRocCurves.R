print("Inside PlottingRocCurves.R")


data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingAllOutputDir()
}else{
  data_type_dir=getTestingAllOutputDir()
}


for(i in 1:length(output_file_names_short))
{
  
  #Plotting ROC curves for Diff. Locations
  current_dataset_name=output_file_names_short[i]
  current_predictions_file=paste(data_type_dir,"/Text_roc_curve_diff_loc_predictions_",current_dataset_name,".txt",sep="")
  predictions_table=t(read.table(current_predictions_file,sep="\t",fill=TRUE))
  rownames(predictions_table)=NULL
  predictions_table=predictions_table[-1,]
  predictions_table_matrix=as.matrix(predictions_table)
  predictions_table_matrix_numeric=predictions_table_matrix
  class(predictions_table_matrix_numeric)='numeric'

  
  current_labels_file=paste(data_type_dir,"/Text_roc_curve_diff_loc_labels_",current_dataset_name,".txt",sep="")
  labels_table=t(read.table(current_labels_file,sep="\t",fill=TRUE))
  rownames(labels_table)=NULL
  labels_table=labels_table[-1,]
  labels_table_matrix=as.matrix(labels_table)
  labels_table_matrix_numeric=labels_table_matrix
  class(labels_table_matrix_numeric)='numeric'
  
  pred.obj=prediction(predictions_table_matrix_numeric,labels_table_matrix_numeric)
  perf.obj=performance(pred.obj,"tpr","fpr")
  
  perf.obj.auc=performance(pred.obj,"auc")
  avg.auc=round(mean(unlist(slot(perf.obj.auc,"y.values"))),2)
  
  fname=paste(data_type_dir,"/Graph_roc_curve_diff_locations_",current_dataset_name,".jpeg",sep="")
  
#   print(fname)
#   inf=slot(perf.obj,"alpha.values")
#   inf_vec=unlist(inf)
#   print("inf_vec")
#   print(head(inf_vec))
  
  
#   if(length(which(predictions_table_matrix_numeric==1))==0)
#   {
#     print("INFINIIIIIIIITEEEE!!! Plot only diagonal")
#     jpeg(fname)
#     plot(perf.obj,main=paste("ROC Curve of Differentiable Locations ",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
#     abline(a=0, b=1,lty=2)
#     dev.off()
#   }else if(length(which(is.infinite(inf_vec)))>0){
#     print("inside else if(length(which(is.infinite(inf_vec)))>0)")
#     jpeg(fname)
#     plot(perf.obj,main=paste("ROC Curve of Differentiable Locations ",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
#     abline(a=0, b=1,lty=2)
#     dev.off()
#   }else{
    print("Regular plot")
    jpeg(fname)
#     plot(perf.obj,avg="threshold",main=paste("ROC Curve of Differentiable Locations Evalucation \n",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
    plot(perf.obj,avg="threshold",main=paste("ROC Curve of Differentiable Locations Evaluation \n"
                                             ,data_type," data for ",current_dataset_name,"\n"
                                             ,"AUC=",avg.auc,sep=""),colorize=FALSE)
    abline(a=0, b=1,lty=2)
    dev.off()
    print("after plot")
#   }
  
  
  
  
  #Plotting ROC curves for Complete Motifs
  current_dataset_name=output_file_names_short[i]
  current_predictions_file_comp_motif=paste(data_type_dir,"/Text_roc_curve_complete_motifs_predictions_",current_dataset_name,".txt",sep="")
  predictions_table_comp_motif=t(read.table(current_predictions_file_comp_motif,sep="\t",fill=TRUE))
  rownames(predictions_table_comp_motif)=NULL
  predictions_table_comp_motif=predictions_table_comp_motif[-1,]
  predictions_table_comp_motif_matrix=as.matrix(predictions_table_comp_motif)
  predictions_table_comp_motif_matrix_numeric=predictions_table_comp_motif_matrix
  class(predictions_table_comp_motif_matrix_numeric)='numeric'
  
  current_labels_comp_motif_file=paste(data_type_dir,"/Text_roc_curve_complete_motifs_labels_",current_dataset_name,".txt",sep="")
  labels_comp_motif_table=t(read.table(current_labels_comp_motif_file,sep="\t",fill=TRUE))
  rownames(labels_comp_motif_table)=NULL
  labels_comp_motif_table=labels_comp_motif_table[-1,]
  labels_comp_table_matrix=as.matrix(labels_comp_motif_table)
  labels_comp_motif_table_matrix_numeric=labels_comp_table_matrix
  class(labels_comp_motif_table_matrix_numeric)='numeric'


  
  pred_obj_comp_motif=prediction(predictions_table_comp_motif_matrix_numeric,labels_comp_motif_table_matrix_numeric)
  perf_obj_comp_motif=performance(pred_obj_comp_motif,"tpr","fpr")
  fname=paste(data_type_dir,"/Graph_roc_curve_complete_motifs_",current_dataset_name,".jpeg",sep="")

  perf_obj_comp_motif_auc=performance(pred_obj_comp_motif,"auc")
  comp_motif_avg_auc=round(mean(unlist(slot(perf_obj_comp_motif_auc,"y.values"))),2)
  


#   inf_comp=slot(perf_obj_comp_motif,"alpha.values")
#   inf_vec_comp=unlist(inf_comp)
#   print("inf_vec")
#   print(head(inf_vec))

# 
# if(length(which(predictions_table_comp_motif_matrix_numeric==1))==0)
#   {
#     print("INFINIIIIIIIITEEEE!!! Plot only diagonal")
#     jpeg(fname)
#     plot(perf_obj_comp_motif,main=paste("ROC Curve of Complete Motifs ",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
#     abline(a=0, b=1,lty=2)
#     dev.off()
# #     quit()
#   }else if(length(which(is.infinite(inf_vec)))>0){
#     print("inside else if(length(which(is.infinite(inf_vec)))>0)")
#     jpeg(fname)
#     plot(perf_obj_comp_motif,main=paste("ROC Curve of Complete Motifs ",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
#     abline(a=0, b=1,lty=2)
#     dev.off()
#   }else{
    print("Regular plot")
    jpeg(fname)
#     plot(perf_obj_comp_motif,avg="threshold",main=paste("ROC Curve of Complete Motifs ",current_dataset_name,",",data_type,sep=""),colorize=FALSE)
    plot(perf_obj_comp_motif,avg="threshold",main=paste("ROC Curve of Complete Motifs Evaluation \n"
                                             ,data_type," data for ",current_dataset_name,"\n"
                                             ,"AUC=",comp_motif_avg_auc,sep=""),colorize=FALSE)    
    abline(a=0, b=1,lty=2)
    dev.off()
    print("after plot")
#   }
}
