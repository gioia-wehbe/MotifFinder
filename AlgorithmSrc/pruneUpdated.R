#Last updated 31-01-15
# print("Inside pruneUpdated.R")

pruneUpdated<-function (object, ...) 
{
  print("inside pruneUpdated()")
  .local <- function (object, nmin, L, gain, C, keep, drop, 
                      state, delete = TRUE) 
  {
    data <- object@data
    cdata <- object@cdata
    A <- alphabet(object)
    cpal <- cpal(object)
    labels <- stlab(object)
    segmented <- object@segmented
    group <- object@group
    if (!missing(gain) & missing(C)) {            #Not entered
      stop(" [!] please provide a cutoff value")
    }
    else if (missing(gain) & !missing(C)) {        #Not entered
      stop(" [!] missing 'gain' argument")
    }
    if (!missing(keep)) 
      {
      if (has.cdata(object)) {
        c.A <- alphabet(object@cdata)
      }
      else {
        c.A <- object@alphabet
      }
      if (!inherits(keep, "stslist")) 
        {
        keep <- seqdef(keep, alphabet = c.A, nr = "#")
#         print("keep")
#         print(keep)
      }
      keep.sl <- seqlength(keep)
#       print("keep.sl")
#       print(keep.sl)
    }
    object <- as(object, "list")
    cnodes <- NULL
    message(" [>] pruning results: ")
    message("   ", format("[L]", width = 5, justify = "right"), 
            format("[nodes]", width = 9, justify = "right"), 
            format("[pruned]", width = 10, justify = "right"))
    
    if (!missing(gain) && is.character(gain)) {       #Not Entered
      if (gain == "G1") {
        gain <- G1
      }
      else if (gain == "G2") {
        gain <- G2
      }
    }
    
    
#     print("length(object)")
#     print(length(object))
    
    for (i in length(object):2) 
      {
#       print("inside for (i in length(object):2) ")
#       print("We are at level --------------------------------------------------------------------------------- i")
#       print(i)
      
      nodes <- object[[i]]
#       print("nodes")
#       print(length(nodes))
      
      parents <- object[[i - 1]]
#       print("parents")
#       print(length(parents))
      
      nbnodes <- unlist(lapply(nodes, function(x) {
        sum(!x@pruned)
      }))   #count the number of nodes that are not pruned
      
#       print("nbnodes")
#       print(nbnodes)
      
      if (!missing(L) && i > (L + 1)) {       #not entered
        nodes <- lapply(nodes, node.prune)
      }
      else 
        {
        if (!missing(keep)) 
          {
#           print("inside if (!missing(keep)) ")
#           print("(i - 1)")
#           print((i - 1))
#           print("max(keep.sl)")
#           print(max(keep.sl))
          if ((i - 1) > max(keep.sl))     #if previous level is >> than the longest possible length of the kept nodes => remove nodes at current level
            {
#             print("inside if ((i - 1) > max(keep.sl)) which prunes")
            nodes <- lapply(nodes, node.prune)
          }
          else 
            {
#               print("inside else which deals with nodes to keep")
            keep.tmp <- keep[keep.sl == i - 1, , drop = FALSE]
            keep.list <- seqconc(keep.tmp)
#             print("cnodes")
#             print(cnodes)
            nodes <- lapply(nodes, nodeKeepModified, keep.list, 
                            clist = cnodes)   #GIOIA: node.keep had a bug too!!!
#             print("after nodeKeepModified")
            
            
          }
        }
        if (!missing(state))    #Not entered
          {
          state.tmp <- seqdecomp(names(nodes))
          state.tmp <- which(rowSums(state.tmp == state) > 
                               0)
          nodes[state.tmp] <- lapply(nodes[state.tmp], 
                                     node.prune)
        }
        if (!missing(nmin))     #Not entered
          {
          nodes <- lapply(nodes, node.nmin, nmin)
        }
        if (!missing(gain))     #Not entered
          {
          nodes <- lapply(nodes, node.gain, plist = parents, 
                          gain = gain, C = C, clist = cnodes)
        }
      }
      
      pruned <- unlist(lapply(nodes, function(x) {
        sum(x@pruned)
      }))
      
#       print("pruned")
#       print(pruned)
      
      plabel <- if (segmented) {
        " node segment(s) pruned.................................................................."
      }
      message("   ", format(i - 1, width = 5), format(sum(nbnodes), 
                                                      width = 9), format(sum(pruned), width = 10))
      
      
      
      
      
if (sum(pruned) > 0) 
        {
#         print("inside if (sum(pruned) > 0) ")
        cnodes <- lapply(nodes, delete.pruned)
#         print("cnodes")
#         print(length(cnodes))
        pruned.id <- which(unlist(lapply(cnodes, function(x) {
          nrow(x@prob) == 0
        })))
#         print("pruned.id")
#         print(pruned.id)
        
        if (length(pruned.id) > 0) 
          {
#           print("inside if (length(pruned.id) > 0) ")
          cnodes <- cnodes[-pruned.id]
        }
        if (delete) 
          {
#           print("inside if (delete) ")
          nodes <- cnodes
#           print("nodes")
#           print(length(nodes))
          if (length(nodes) > 0) 
            {
#             print("inside if (length(nodes) > 0) ")
            remaining <- lapply(nodes, function(x) {
              rownames(x@prob)
            })
            rplist <- unlist(lapply(nodes, node.parent))
            parents <- lapply(parents, set.leaves, remaining, 
                              rplist)
          }
          else {
            parents <- lapply(parents, function(x) {
              x@leaf[] <- TRUE
              x
            })
          }
        }
      }
      else if (!delete) {
#         print("inside else if (!delete) ")
        cnodes <- nodes
      }
      else if(delete)#GIOIA: Bug in original function: did not consider this case!!!
      {
#         print("inside else if(delete)#GIOIA: Bug in original function: did not consider this case!!!")
        cnodes <- nodes
      }
      else
      {
#         print("EEEEEEEEEEEEEEEEEEEEEEElllllllseeee... did not initialize cnodes!!")
#         print("nodes")
#         print(length(nodes))
      }






      if (length(nodes) == 0) 
        {
#         print("inside if (length(nodes) == 0)")
        object <- object[-i]
      }
      else {
#         print("inside else if (length(nodes) == 0) ")
        object[[i]] <- nodes
      }
      object[[i - 1]] <- parents
    }#end of loop: for (i in length(object):2) 
    
#     print("final pruned object:")
#     print(object)

    object <- new("PSTf", object, data = data, cdata = cdata, 
                  alphabet = A, cpal = cpal, labels = labels, segmented = segmented, 
                  group = group, call = match.call(), logLik = as.numeric(NULL))
    debut.lik <- Sys.time()
    return(object)


  }
  .local(object, ...)
}


environment(pruneUpdated)=asNamespace("PST")



