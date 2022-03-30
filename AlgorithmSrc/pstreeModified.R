##################################################################################################
# A function that builds the probabilistic suffix trees based on relative frequency per position
# Source: PST library
# Modified by: Gioia Wehbe
# Last updated: 20-03-15
##################################################################################################

pstreeModified<-function (object, group, L, ...) 
{
  print("inside pstreeModified")
  
  .local <- function (object, group, L, cdata = NULL, stationary = TRUE, 
                      nmin = 1, ymin = NULL, pmin=0, weighted = TRUE, with.missing = FALSE) 
  {
    
    
    
    #************************************************************************
    # Initializaion of needed variables and validating parameters....
    #************************************************************************
    
#     library("plyr")
    debut <- Sys.time()
    if (!stationary & !flist("pstree", "stationary")) {   #Not needed
      stop(" [!] 'stationary=FALSE' not implemented", call. = FALSE)
    }
    if (!is.null(cdata) & !flist("pstree", "cdata")) {    #Not needed
      stop(" [!] argument cdata not allowed", call. = FALSE)
    }
    if (missing(L)) {     #This means that longest possible motif to consider is seqlength-1. That's why we add $ to indicate end of string
      L <- max(seqlength(object)) - 1
    }
    if (!weighted || is.null(attr(object, "weights"))) {    #not needed
      attr(object, "weights") <- rep(1, nrow(object))
    }
    if (all(attr(object, "weights") == 1)) {    #not needed
      weighted <- FALSE
    }
    A <- alphabet(object)   #possible SNP values: 0,1,2
    StCol <- cpal(object)   #state color: color that corresponds to each SNP value (NOT NEEDED)
    StLab <- stlab(object)  #state label: a label for each SNP (0,1 and 2 on our case)
    GroupCol=NULL    #GIOIA: Initialization of list of colors per group
    GroupLab=NULL    #GIOIA: Initialization of list of labels per group
    sl <- seqlength(object)
    if (with.missing)       #Not needed
    {   #incase some of the states are missing (not the case...)
      A <- c(A, attr(object, "nr"))
      StCol <- c(StCol, attr(object, "missing.color"))
      StLab <- c(StLab, "missing")
    }
    names(StCol) <- A
    
    if (!missing(group)) 
      {
      GroupLab <- levels(factor(group))   #GIOIA: create list of labels indicating the different groups

      num_labels=length(GroupLab)

#       GroupCol <- palette(rainbow(num_labels,start=0.1,alpha=1))  #GIOIA: create list of colors based on the number of different groups
      rainbow=rainbow(num_labels,start=0.1,alpha=1)      

      GroupCol <- rainbow

      
      names(GroupCol)=GroupLab      #GIOIA: set the index of each color in GroupCol to Group label
      
      if (!is.factor(group)) 
        {
        group <- factor(group)
      }
      nbgroup <- length(levels(group))
      segmented <- TRUE     #Not Needed
      message(" [>] ", nbgroup, " groups")
      }
    else {      #Not needed
      group <- factor(NULL)
      segmented <- FALSE
    }
    
    nodes.list <- vector("list", length = L + 1)    # the list that will include all nodes in the form of a tree
    #where every index includes the list of node of the length of the index
    
    message(" [>] ", nrow(object), " sequence(s) - min/max length: ", 
            min(sl), "/", max(sl))
    message(" [>] max. depth L=", L, ", nmin=", nmin, if (!is.null(ymin)) {
      paste(", ymin=", ymin, sep = "")
    })
    message("   ", format("[L]", width = 5, justify = "right"), 
            format("[nodes]", width = 9, justify = "right"))
    
    
    #****************************************************************************************************
    # start of the algorithms
    #****************************************************************************************************
    for (i in 0:L) #for every possible length of a motif starting from 1 till the whole motif
      {

      

      if (segmented)  #if there are multiple classes (Always the case)
        {
        print("if (segmented)")
        tmp <- NULL       #initialize list of all highly probable motifs per position of all possible lengths
        
        for (g in 1:nbgroup)        #for every possible group the motif can belong to
          {
          
#           print("for (g in 1:nbgroup) ")
#           print("current group")
#           print(g)
          data <- object[group == levels(group)[g], ]     #get the data
          
          tmp.cdata <- if (!is.null(cdata)) {         
            cdata[group == levels(group)[1], ]
          }     
          else {
            NULL
          }                                  
          

          #get the frequency per position of every possible motif of length L for the current class g
          #however!! add only those motifs with probability per position>>pmin (GIOIA: parameter added by me)
#           print("data")
#           print(data)
#           print("current L")
#           print(i)
          ccounts <- suppressMessages(cprobModified(data, L = i, 
                                            cdata = tmp.cdata, stationary = stationary, 
                                            nmin = nmin, pmin=pmin, prob = FALSE, weighted = weighted, 
                                            with.missing = with.missing, to.list = TRUE))
#           print("after cprobModified")

          if(!is.null(ccounts))
          {
#             print("inside if(!is.null(ccounts))")
            
            #GIOIA: add the positions at which the motif of length l was found to ccount
            ccounts <- lapply(ccounts, function(x)
            {
              pos_start_ind=length(A)+3
              pos_end_ind=ncol(x)
              positions=t(as.matrix(x[,pos_start_ind:pos_end_ind]))
              cnames=colnames(positions)
              if(length(cnames)==0)
                cnames=1
              cnames=paste("p.",cnames,sep="")
              colnames(positions)=cnames
              x=t(as.matrix(x[,1:pos_start_ind-1]))
              rownames(x)=NA
              cbind(x,group = g, position = as.integer(rownames(x)),positions)
            } )
            
#             print("after  ccounts <- lapply(ccounts, function(x)")
            
            tmp <- merge.cprob.modified(tmp, ccounts)   #add counts of motifs of length L and group g to the list of 
                                                        #all motifs of length L repeat for the next possible group g
#             print("after merge.cprob.modified")
          }
          else
          {
#              print("no motifs with high relative freq per position were foun for current length and current group")
          }
#######################################################################################3          
        }#end of group loop                           
      }
      else      #if not segmented: never the case since there are always classes (Not needed)
        {
        tmp <- suppressMessages(cprob(object, L = i, 
                                      cdata = cdata, stationary = stationary, nmin = nmin, 
                                      prob = FALSE, weighted = weighted, with.missing = with.missing, 
                                      to.list = TRUE))
        tmp <- lapply(tmp, function(x) {
          cbind(x, group = NA, position = as.integer(rownames(x)))
        })
      }
      
      nodes.names <- names(tmp)   #names of all possible motifs of length L that have high frequency per position for at least one of the groups

      if(!is.null(nodes.names))   #GIOIA: If motifs were kept (there are motifs where prob > pmin) => create tree nodes to add to tree
      {
#         print("inside  if(!is.null(nodes.names))  ")
        tmp.list <- lapply(seq_len(length(tmp)), function(n)
          { # convert every such motif to a PSTr node which is a probabilistic suffix tree node that represents 
          # the motif and the details about it including its positions, its max frequency per position, the probability
          # of it being followed by each of the states....
          
          current_ccount=tmp[[n]]
          index_start_ind=length(A)+2
          index_end_ind=ncol(current_ccount)
          curr_index_matrix=tmp[[n]][,index_start_ind:index_end_ind]
          if(class(curr_index_matrix)!="matrix")
                curr_index_matrix=t(as.matrix(curr_index_matrix))
          index_cnames=colnames(curr_index_matrix)
          new("PSTr", 
              path = nodes.names[n],    #The motif
              counts = tmp[[n]][, A, drop = FALSE],     #The counts to which the probability distributions of each of the states following the motif
              n = tmp[[n]][, "n", drop = FALSE],     #The number of occurrences of the context (motifs) in the learning sample 
              order = i,                              #The depth of the node in the tree which is equivalent to the length of the motif
              ymin=ymin,                              #Not needed
              index = tmp[[n]][, index_cnames, drop = FALSE])     #it includes other details about the node including its positions
        } )
        
#         print("after tmp.list <- lapply(seq_len(length(tmp)), function(n)")
        
        nbnodes <- length(tmp.list)   #number of nodes
        names(tmp.list) <- nodes.names    #set the names of the PSTr node list to the names of the motif each node  
                                          # corresponds to
        if (i > 0) #for all nodes except the starting node 'e'
          {
#           print("inside if (i > 0)")
          parents <- nodes.list[[i]] #create a list 'parents' which includes the parents of the current length-i-nodes 
                                      #(temp.list). parent nodes are of length i-1
          
          child.list <- lapply(tmp.list, function(x) {
            rownames(x@prob)
          })      #create a list 'child.list' which is the row.names of the nodes in temp.list (which are contexts of lengths i+1)
          
          #create a list 'rplist' which includes the names of the parents of nodes in temp.list (which are contexts of length i-1) 
          rplist <- unlist(lapply(tmp.list, node.parent))   #names of the nodes in the list 'parents'
          #
          nodes.list[[i]] <- lapply(parents, set.leaves, 
                                    child.list, rplist)     #set the leaves of the list 'parents' to the new child node
#         print("after set.leaves")
          }   #of length i found (child.list). and reinitialize node.list[i] to the new parents nodes (of length i-1) with their leaves set
        message("   ", format(i, width = 5), format(nbnodes, 
                                                    width = 9))

        nodes.list[[i + 1]] <- tmp.list     #Add the new nodes found (of length i) to the node.list (at index i+1 ) since i starts at 0
        
      }
      else    #GIOIA: If no motifs were kept (no motifs where prob> pmin) => This is the longest motif with high relative frequency per position
      {# So end tree
#         print("no differentiable motifs found at this level...")
                break;
      }
    }
    nodes.list[[i]] <- lapply(nodes.list[[i]], node.leaf)   #GIOIA: set the longest nodes last found as leaves in the tree to tree
    
    #Remove extra empty locations from nodes.list
#     start_indx=i+1
#     print("start_indx")
#     print(start_indx)
#     end_indx=L+1
#     print("end_indx")
#     print(end_indx)
#     print("nodes.list")
#     print(nodes.list)
    nodes.list=nodes.list[!sapply(nodes.list,is.null)]
#     print("nodes.list after eliminating null")
#     print(nodes.list)
#     nodes.list=nodes.list[-start_indx:-end_indx]
#       print("nodes.list")
#   print(nodes.list)
    
    if (is.null(cdata)) #Not needed
      {
      cdata <- object[-(1:nrow(object)), ]
    }
    else {               #Not needed
      if (with.missing) {
        cdata <- seqdef(cdata, nr = "#", alphabet = c(alphabet(cdata), 
                                                      attr(cdata, "nr")), labels = c(stlab(cdata), 
                                                                                     "missing"), cpal = c(cpal(cdata), attr(cdata, 
                                                                                                                            "missing.color")), xtstep = attr(cdata, "xtstep"))
      }
    }
    
    # Build the final tree which is a PSTf object
#     The class "PSTf" is the flat representation of a probabilistic suffix tree (PST) storing a variable
#     length Markov chain model. The flat representation is a list where each element corresponds to
#     a given depth. It is the prefered representation and is used by all functions for model fitting and
#     sequence analysis with PST. The nested representation "PSTr" is used only for printing and plotting
#     PSTs.
#     
    res <- new("PSTf", nodes.list, data = object, cdata = cdata, 
               alphabet = A, cpal = GroupCol, labels = GroupLab, segmented = segmented, 
               group = group, call = match.call(), logLik = as.numeric(NULL))
    #GIOIA: eliminated calculation of likelihood since not needed and time consuming
    return(res)
  }
  .local(object, group, L, ...)
}

environment(pstreeModified)=asNamespace("PST")
