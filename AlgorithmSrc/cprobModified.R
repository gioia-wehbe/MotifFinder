##################################################################################################
# A function that computes the relative frequencies per position for motifs of a specific length L
# in each group. It returns those motifs whose relative frequency per position>pmin
# Source: PST library
# Modified by: Gioia Wehbe
# Last updated: 20-03-15
##################################################################################################

cprobModified<-function (object, L, ...) 
{
  print("inside cprobModified")
  .local <- function (object, L, cdata = NULL, context, stationary = TRUE, 
                      nmin = 1,  pmin=0,prob = TRUE, weighted = TRUE, with.missing = FALSE, 
                      to.list = FALSE) 
  {
#     library("plyr")
    debut <- Sys.time()   #for timing 
    if (!is.null(cdata) & !flist("cprob", "cdata"))     #Not needed
      {
      stop(" [!] argument cdata not available", call. = FALSE)
    }
    statl <- alphabet(object)   #List of states ie SNP values (0, 1, 2)
    
    nr <- attr(object, "nr")    #Not needed
    if (with.missing)           #Not needed
      {
      statl <- c(statl, nr)
    }
    
    sl <- seqlength(object)     #sequence length
    sl.max <- max(sl)           #Max sequence length (in our case, all the same length)
    
    if (L > (sl.max - 1))       # validating that the length of the motifs under search is valid (within the length of the sequences)
      {
      stop(" [!] sequence length <= L")
    }
    
    nbseq <- nrow(object)   #number of sequences in the data
    weights <- attr(object, "weights")      #Not needed
    if (!weighted || is.null(weights))        #Not needed
      {
      weights <- rep(1, nrow(object))
    }
    object <- as.matrix(object)            #convert data to matrix

    snp_locations=colnames(object)      #Gioia
    if (!missing(context))          #validating that the motif under search is a valid motif NOT NEEDED since we are not specifying motif... we are only specifying motif length
      {
      tmp <- seqdecomp(context)
      if (any(!tmp %in% statl) & context != "e") 
        {
        stop(" [!] one or more symbol in context not in alphabet")
      }
      L <- ncol(tmp)        #set L to the length of the context in search
    }
    
    
    message(" [>] ", nbseq, " sequences, min/max length: ", 
            min(sl), "/", max(sl))            #displays the number of sequences
    
    states <- factor(object[, (L + 1):sl.max], levels = statl)    #all possible states that might follow a context of length L in the current data 
    states_matrix <- as.matrix(object[, (L + 1):sl.max])        #states in a matrix

    if(nbseq==1)
      states_matrix=t(states_matrix)
    

    
    contexts <- matrix(nrow = nbseq, ncol = sl.max - L)   #declaration of a matrix of all possible contexts of length L found in the give data

    #Initialization of the contexts matrix by shifting the data 
    if (L == 0)     
      {
      contexts[] <- "e"
    }
    
    
    else 
      {
      if (is.null(cdata)) #Not needed
        {
        cdata <- object
      }
      else        #not needed
        {
        cdata <- as.matrix(cdata)
      }
      
      
      for (p in (L + 1):sl.max) 
      {
        contexts[, p - L] <- cdata[, (p - L)]
        if (L > 1) 
        {
          for (c in (L - 1):1) 
          {
            contexts[, p - L] <- paste(contexts[, p -L], cdata[, (p - c)], sep = "-")
          }
        }
      }
    }#end of initialization of the contexts matrix
    
    contexts_matrix=contexts    #GIOIA: contexts in the form of a matrix (every column is a list of all possible contexts at certain position in the current data)

    colnames(contexts_matrix)=snp_locations[1:ncol(contexts_matrix)]
    contexts <- as.vector(contexts)   #all contexts stored in a vector/list (used to calculate over all frequency of a context irrespective of the position)    #From old code not very much needed in our code: Just in case missing was set to TRUE. NOT NEEDED

    
    if (!missing(context)) 
    {
      sel <- contexts == context
      contexts <- contexts[sel]
      states <- states[sel]
      weights <- weights[sel]
    }
    
    message(" [>] computing prob., L=", L, ", ", length(unique(contexts)), 
            " distinct context(s)")
    if (stationary)   #Always the case.... always stationary
    {
      #initialization of the final result matrix
      res=matrix(ncol=length(statl)+3)    
      colnames(res) <- c(statl, "n","rel_n","pos")
      positions=matrix(ncol=1)      #Gioia: a matrix that includes the positions of the contexts
      colnames(positions)=c("pos")
     
      #n is a list of the total counts of all possible contexts of length L (followed by states) in the data irrespective of their position
      n <- rowSums(xtabs(~contexts + states)[, , drop = FALSE])
      #Start the filtration per position
      #GIOIA: pstn is the position which is the number of columns in the contexts matrix
      actual_positions=colnames(contexts_matrix)
      

      
      for(pstn in 1:dim(contexts_matrix)[2]) 
      {
        #get all contexts at position pstn and convert to vector (so we can use with xtabs)
        current_actual_position=actual_positions[pstn]
        current_actual_position=as.numeric(current_actual_position)

        curr_contexts=as.vector(contexts_matrix[,pstn]) 
        #get all the states that might follow contexts at position pstn (ie get all states that are found at position 
        #pstn+1 in the data set)

        curr_states=states_matrix[,pstn]
        #factor current states so we can use with xtabs()
        curr_states=factor(curr_states,levels=statl)
        #freq per position:  compute the frequency of occurence of every context at position pstn followed by the 
        #possible states (at position pstn+1) ONLY IN THE CURRENT POSITION!!! NOT IN ALL DATA!!!

        freq <- xtabs(weights ~ curr_contexts + curr_states)[, , drop = FALSE]
        #the relative frequency normalized by the number of sequences (for the current group) so that we can compare 
        #Gioia: groups irrespective of the different number of sequences in each one 
        rel_freq=freq/nbseq       
        
        if (prob)#Not needed in our case since we don't want to output probability (Prob=FALSE)
        {
          freq <- freq/rowSums(freq) 
        }
        
        #Gioia: n_pos is a list of the total counts of all possible contexts of length L in the data AT THE CURRENT POSITION!!!
        n_pos <- rowSums(xtabs(~curr_contexts + curr_states)[, , drop = FALSE])
        #GIOIA: it is the relative frequency of contexts of length L normalized by the number of sequences (in the current 
        #group) so that we can compare between groups with different numbers of sequences
        rel_n=n_pos/nbseq

        pstns=rep(current_actual_position,length(freq))  #GIOIA: a list of the position of every context
        curr_res <- cbind(freq, n_pos,rel_n,pstns)  #GIOIA: result that stores all frequencies PER POSITION!! 
        #GIOIA: result that stores RELATIVE FRQUENCIES (normalized by the number of sequences in the current group data)
        rel_res=cbind(rel_freq,rel_n,pstns)         
        #GIOIA: Compares the relative frequency (normalized by the number of sequences in the current group data) 
        #per position to the minimum frequency specified in the parameters (TODO: Try various pmin parameters)
        #Gioia:pmin was added in the parameters(default = 0 and thus doesn't influence anything if decided not to use)
        if (pmin > 0)     
        {
          #Gioia: find which contexts of length L have a relative frequency per position less than pmin
          pmin.del <- which(rel_res[, "rel_n"] < pmin)
          if (length(pmin.del) > 0)      #Gioia: if such contexts were found
            {
            #Gioia: remove these context from the count result of the current position
            curr_res <- curr_res[-pmin.del, , drop = FALSE]   
            #Gioia: output how many contexts were removed...
            message(" [>] removing ", length(pmin.del),     
                      " context(s) where proba<", pmin)
            }
        }#Gioia: end of pmin >0
        
        # GIOIA: add the contexts with relative freq per position > pmin to the final result
        if(nrow(curr_res)!=0)    #GIOIA: if such context exist in the current position:
        {
          if(is.na(res[1,1])) #GIOIA: if the final result is still empty
          {
            #Gioia: set the first row of the final result to the first row of the current result that has the high 
            #frequency context (that is to replace the NAs that are originally found in res)
            res[1,]=curr_res[1,]  
            #GIOIA: if the results of the current positon are more than one (which has replaced the NAs)
            if(nrow(curr_res)>1)    
            {
              #GIOIA: bind all the rows of the current result to the final result except for the first one since it has 
              #already been added to replace the NAs
              res=rbind(res,curr_res[-1,])  
            }
            #GIOIA: name the rows of the final result same as the current result (copy the names of the contexts to 
            #the final result)
            rownames(res)=rownames(curr_res)  
          }#End of 'If the final result is still empty'
          else  #in case this is not the first time we add to the final result(in case the final result is not empty)
            {
              res=rbind(res,curr_res)
            }
          }# End of 'if highly frequent contexts exist in the current position'
         else   #GIOIA: in case the current position doesn't have any highly frequent motifs (curr_res is empty) 
         { #don't add to the final result at all
#             print("curr_res is empty dont add")
         }
      }#End of for loop. looked at all positions and collected all highly frequent motifs


      #GIOIA: the same context might appeare more than once since its found with high frequency in different positions
      #we eliminate redundancies by selecting the context with maximum relative Frequency (check maxRelN) 
      if(!is.na(res[1,1]))   # GIOIA: if the final result is not completely empty
      {
        #selects the context with the highest relative frequency to keep as a node, but keeps all positions
        res=sapply(by(res,rownames(res),maxRelN),identity)
        res=t(res)    #transpose the result matrix since the previous function turns rows into columns...
        
        #Finally remove contexts with counts<nmin in case nmin is specified  (never the case) #not needed
        if (nmin > 1) {
          nmin.del <- which(res[, "n"] < nmin)
          if (length(nmin.del) > 0) {
            res <- res[-nmin.del, , drop = FALSE]
            message(" [>] removing ", length(nmin.del), 
                    " context(s) where n<", nmin)
          }
        }
      }
      else
      {
#         print("result is completely empty")
        return(NULL)
      }
      #Gioia: convert all possible positions of each node from a comma-dellimited charater string to a numeric matrix
      positions=as.matrix(res[,'pos'])    #Gioia: get the positions as a list of strings from res
      #Gioia: assign the name of the context that corresponds to each set of positions 
      rownames(positions)=rownames(res)   
      colnames(positions)="pos" # name the column of the position matrix
      positions_list=as.list(positions[,1])   #Gioia: get the list version of the position matrix
      #Gioia split the comma-delimited character strings to vectors and we get  a list of vectors
      position_list_split=sapply(positions_list,function(x){strsplit(as.character(x),',',fixed=TRUE)})
      #Gioia: convert it to a list of one-row matrices instead of vectors to be able to use with rbind.fill.matrix
      position_list_split=lapply(position_list_split,function(x){t(as.matrix(x,nrow=1,ncol=length(x)))})
      #Gioia: convert into a numeric position matrix where each element indicates a position and each row is a 
      # context. NB: if one context has less positions compared to the others, the empty locations are filled with NA
      pos_matrix=rbind.fill.matrix(position_list_split)
      #Gioia: fix the rownames of the positions matrix
      rownames(pos_matrix)=rownames(positions)
      #Gioia: set mode to numeric
      mode(pos_matrix)='numeric'
      
      #Gioia: add the position matrix of the end of the result matrix making sure that the mode is numeric
      c_names=colnames(res)[1:ncol(res)-1]
      res=matrix(res[,1:ncol(res)-1],nrow=nrow(res),dimnames=list(rownames(res),c_names))
      res=as.matrix(res)
      mode(res)='numeric'
      res=cbind(res,pos_matrix)
    } #End of if stationary
    else #if not stationary: this is never the case. NOT NEEED
      {
      t <- (L + 1):sl.max
      pos <- matrix(t, ncol = length(t), nrow = nbseq, 
                    byrow = T)
      pos <- as.vector(pos)
      pos <- factor(pos)
      if (!missing(context)) {
        pos <- pos[sel]
      }
      tmat <- xtabs(weights ~ pos + states + contexts)
      n <- xtabs(~pos + states + contexts)
      context.list <- dimnames(tmat)$contexts
      res <- lapply(context.list, function(idx) {
        if (nrow(tmat[, , idx, drop = FALSE]) == 1) {
          freq <- t(as.matrix(tmat[, , idx]))
          rownames(freq) <- rownames(tmat[, , idx, drop = FALSE])
          fs <- sum(n[, , idx])
        }
        else {
          freq <- tmat[, , idx]
          fs <- rowSums(n[, , idx])
        }
        if (prob) {
          freq <- freq/rowSums(freq)
        }
        pplist <- cbind(freq, n = fs)
        pplist <- pplist[order(as.numeric(rownames(pplist))), 
                         , drop = FALSE]
        nmin.del <- which(pplist[, "n"] < nmin)
        if (length(nmin.del) > 0) {
          pplist <- pplist[-nmin.del, , drop = FALSE]
        }
        return(pplist)
      })
      names(res) <- context.list
      if (L > 0) {
        empty <- which(unlist(lapply(res, function(x) {
          is.null(x) || nrow(x) == 0
        })))
        if (length(empty) > 0) {
          res <- res[-empty]
        }
      }
    }
    
    

    if (L > 0 & !with.missing) #does not affect us... NOT NEEDED
      {
      if (nr %in% c("?", "*")) {
        nr <- paste("\\", nr, sep = "")
      }
      hasMiss <- if (stationary) 
        {
        grep(nr, rownames(res))
      }
      else {
        grep(nr, names(res))
      }
      if (length(hasMiss) > 0) {
        message(" [>] removing ", length(hasMiss), " context(s) containing missing values")
        res <- if (stationary) {
          res[-hasMiss, , drop = FALSE]
        }
        else {
          res[-hasMiss]
        }
      }
    }
    
    
    fin <- Sys.time() #output time
    message(" [>] total time: ", format(round(fin - debut, 3)))
    
    #convert the res result from a matrix to a list where each index is a context's details and is indexed by the 
    #context's name. this is used if to.list is TRUE. This is always the case when used to build the pstree
    if (stationary & to.list) 
      {
      res <- lapply(1:nrow(res), function(i) res[i, , drop = FALSE])
      nodes.names <- unlist(lapply(res, rownames))
      names(res) <- nodes.names
      res <- lapply(res, function(x) {rownames(x) <- NA 
                                      x})
    }
    return(res)
  }
  .local(object, L, ...)
}

environment(cprobModified)=asNamespace("PST")
