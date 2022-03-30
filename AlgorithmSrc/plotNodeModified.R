#########################################################################
# Plots nodes
# calls: plotNodeProbModified
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 31-01-15
#########################################################################



#Source: https://r-forge.r-project.org/scm/viewvc.php/pkg/R/plotNode.R?view=markup&root=pst&pathrev=173 
plotNodeModified <- function(x1, x2, subtree, seglist, nPar, horiz = FALSE, gratio, max.level, cex, nc, cpal) 
{
  scale <- seq(0, 1, 0.2)
  
  ## Retrieving requested attributes

  prob <- subtree@prob

  index <- subtree@index
  nbseg <- nrow(seglist)
  seg.lab <- rownames(seglist)
  has.segment <- rownames(seglist) %in% rownames(index)  
  
  max_rel_n=0
  majority_group=NULL
  for (g in 1:nbseg) #for every group
  {
    ##########################################Gioia
    if (has.segment[g])   #if the current node belongs to this group...
    {
      idseg <- rownames(seglist)[g]
      rel_n=index[idseg,"rel_n"]
      if(rel_n>max_rel_n)
      {
        majority_group=g
        max_rel_n=rel_n
      }
    }
  }
  majority_group_color=cpal[majority_group]
  path <- subtree@path
  n<-subtree@n
  state <- seqdecomp(path)[1]
  alphabet <- colnames(prob)
  depth <- subtree@order
  pruned <- subtree@pruned
  if (is.null(pruned)) {pruned <- FALSE}
  node.id <- seqdecomp(path)[1]
  children <- which.child(subtree)
  inner <- !length(children)==0 && x1 != x2 && !depth==max.level
  
  bx <- plotNodeLimit(x1, x2, subtree, max.level)
  X <- bx$x
  Y <- subtree@order
  
  ## Setting node parameters
  i <- if (inner) 1 else 2
  
  node.type <- Xtract("node.type", nPar, default = c("prob", "prob"), i)
  node.size <- Xtract("node.size", nPar, default = 0.6)

  Node.lim <- ((node.size/2)*gratio)
  pfactor <- Xtract("pfactor", nPar, default = 0)
  node.size <- node.size*(1+pfactor)^(max.level-depth)
  ns.adj <- ((node.size/2)*gratio)
  
  ## color palette for the context id
  c.cpal <- Xtract("c.cpal", nPar, default = cpal)
  pch <- Xtract("pch", nPar, default = 1L:2, i)
  pruned.col <- Xtract("pruned.col", nPar, default = "red")
  root.col <- Xtract("root.col", nPar, default = "grey")
  bg <- Xtract("bg", nPar, default = par("bg"), i)
  pform <- Xtract("pform", nPar, default = "SPS")
  fg.col <- Xtract("fg.col", nPar, default = "grey")
  fg.del.col <- Xtract("fg.del.col", nPar, default=fg.col)
  
  lab.col <- Xtract("lab.col", nPar, default = par("col"), i)
  lab.cex <- Xtract("lab.cex", nPar, default = cex)
  lab.font <- Xtract("lab.font", nPar, default = par("font"), i)
  lab.type <- Xtract("lab.type", nPar, default = NULL, i)
  lab.srt <- Xtract("lab.srt", nPar, default = 0, i)
  lab.pos <- Xtract("lab.pos", nPar, default = NULL, i)
  lab.offset <- Xtract("lab.offset", nPar, default = 0.5, i)
  c.size <- Xtract("c.size", nPar, default = (node.size/2)*0.66*gratio)
  
  t.col <- Xtract("t.col", nPar, default = "black", i)
  t.cex <- Xtract("t.cex", nPar, default = lab.cex)
  t.font <- Xtract("t.font", nPar, default = par("font"), i)
  
  leave.csize <- Xtract("leave.csize", nPar, default=c.size*0.5)
  leave.lh <- Xtract("leave.lh", nPar, default=0.1)
  leave.lw <- Xtract("leave.lw", nPar, default=node.size)
  lty <- Xtract("lty", nPar, default = par("lty"), i)
  lwd <- Xtract("lwd", nPar, default = 3, i)
  
  stem.height <- Xtract("stem.height", nPar, default=0.4)
  
  ## if verbose==TRUE
  if (getOption("verbose")) {
    cat(" [v] node: ", path, " - L=", Y, if (inner) { " (inner node)" } else { "(leaf)" }, "\n" )
    if (!is.null(nPar)) {
      cat("    [-] node.size", node.size, "\n")
    }
    cat("    [-] x1=", x1, ", x2=", x2, " --> X=", X, "\n")
  }
  
  
  ## TEXT APPEARING in THE NODE Gioia: text appearing in the square, this should be the class proba! We eliminated 
  # the text to keep it simple (we can add if we have small trees) 
  nodeText <- NULL
  
  if (!is.null(lab.type)) 
  {
    nodeText <- vector("character", length=length(lab.type))
    for (ltidx in 1:length(lab.type)) 
    {
      if (lab.type[ltidx]=="path") #don't want path
      {
        if (pform=="SPS" & depth>1) path <- suppressMessages(seqformat(path, from="STS", to="SPS", 
                                        SPS.out=list(xfix="", sdsep="/"), compressed=TRUE))
        nodeText[ltidx] <- path
      } 
      else if (lab.type[ltidx]=="state") #don't want state
      { 
        nodeText[ltidx] <- state 
      }
      else if (lab.type[ltidx]=="n") 
      { 
        nodeText[ltidx] <- round(sum(subtree@n, na.rm=TRUE),1) 
      }
      else if (lab.type[ltidx]=="prob") 
      { 
        nodeText[ltidx] <- paste("(",paste(round(prob,2),collapse=","),")",sep="")
      }
      else if (lab.type[ltidx]=="rel.n") 
      { 
        nodeText[ltidx] <- paste("(",paste(round(index[,"rel_n"],2),collapse=","),")",sep="")
      }
    }
    nodeText <- paste(nodeText, collapse="\n")
  }
  

  
  
  
  if (horiz) #Never the case for now  Not needed
  {
    tmp <- Y
    Y <- X
    X <- tmp
    offset <- nchar(nodeText)/2
    inches <- node.size/2
  } 
  else 
  {
    offset <- nchar(nodeText)/2
    inches <- FALSE
  }
  
  if (node.type=="prob") #Type is alway prob... it is never path
  {
    ## Plotting the circle with the node id   Gioia: Plotting the circle with state in it:
    
    #1- Selecting circle color:
    id.fg <- if (all(pruned)) { fg.del.col } else { fg.col }  #Color if pruned (I think)
    
    #Color if not pruned (based on the color of the state): Gioia: This should be based on the color of the classes to which this context belongs
    #     id.bg <- if (node.id=="e") { root.col } else { c.cpal[which(node.id==names(c.cpal))] }  
    id.bg <- if (node.id=="e") { root.col } else { "grey" }  
    
    if (horiz) #Never the case for now
    {
      ## Middle of the edge stemming from children
      ## yMid <- (X-((node.size/2)*gratio)+X-0.5)/2
      yMid <- (X-((node.size/2)*gratio)+X-(1-stem.height))/2
      ## segments(yMid, Y, X, Y, col=fg.col, lty=lty, lwd=lwd)
      symbols(yMid, Y, circles=c.size, bg=id.bg, fg=id.fg, inches=nc/4, add=TRUE)
      text(yMid, Y, asTxt(node.id), cex=t.cex, font=t.font, col=t.col)
    } 
    else #If vertical
    {#Plots the circle with the label (label of circle should still be state)
      yMid <- (Y-((node.size/2)*gratio)+Y-(1-stem.height))/2
      ## segments(X, yMid, X, Y, col=fg.col, lty=lty, lwd=lwd)
      symbols(X, yMid, circles=c.size, bg=id.bg , fg=id.fg, add=TRUE, inches=FALSE)
      text(X, yMid, asTxt(node.id), cex=t.cex, font=t.font, col=majority_group_color)
    }
    
    Node.ytop <- if (horiz) {Y-(node.size/2)} else { Y-(node.size/2)*gratio }
    Node.ybottom <- if (horiz) {Y+(node.size/2)} else { Y+(node.size/2)*gratio }
    Node.xleft <- if (horiz) { X+((node.size/2)*gratio) } else { X-(node.size/2) }
    Node.xright <- if (horiz) { X-((node.size/2)*gratio) } else { X+(node.size/2) }
    
    if (nrow(prob)>1) 
    {
      if (path=="e") {
        if (horiz) { probAxes <- c("bottom", "right") } else { probAxes <- c("top", "left") }
      } else if (!inner) {
        if (horiz) { probAxes <- c("no", "no") } else { probAxes <- c("bottom", "no") }
      } else { probAxes <- c("no", "no") }
    } 
    else {
      probAxes <- if (path=="e") {
        if (horiz) { c("no", "right") } else { c("no", "left") } } 
      else {c("no", "no")}
    }
    
    ## 
    if (!inner) #if current node is an outer nodes
    {
      ## Bare verticale en dessous du rectangle : this must mean vertical line under the rectangle
      if (horiz) 
      {
        segments(X, Y, X+Node.lim+leave.lh, Y, col=fg.col, lty=lty, lwd=lwd)
      } 
      else 
      {
        segments(X, Y, X, Y+Node.lim+leave.lh, col=fg.col, lty=lty, lwd=lwd)
      }
      
      if (all(subtree@leaf, na.rm=TRUE)) #If current node is a leaf
      {
        ## Leaf indicator : the line that indicates that this outer node is a leaf
        if (horiz) 
        {
          ## Bare horizontale du rateau
          segments(X+Node.lim+leave.lh, Y-(leave.lw/4), 
                   X+Node.lim+leave.lh, Y+(leave.lw/4), 
                   col=fg.col, lty=lty, lwd=lwd)
        } 
        else 
        {
          ## Bare horizontale du rateau
          segments(X-(leave.lw/4), Y+Node.lim+leave.lh, X+(leave.lw/4), 
                   Y+Node.lim+leave.lh, col=fg.col, lty=lty, lwd=lwd)
        }
      }
      else 
      {
        if (horiz) 
        {
          symbols(X+Node.lim+leave.lh, Y, circles=1, inches=(nc/4)*0.5, add=TRUE,
                  fg=fg.col, bg=fg.col)
        } else 
        {##################################################################
          symbols(X, Y+Node.lim+leave.lh, circles=leave.csize, inches=FALSE, add=TRUE,
                  fg=fg.col, bg=fg.col)
        }
      }
    }
    
    plotNodeProbModified(Node.xleft, Node.ybottom, Node.xright, Node.ytop, prob=prob,n=n, seglist=seglist, state=NULL, 
                         cpal=cpal, pruned=pruned, index=index, axes=probAxes)    #Gioia: this is the function that plots the square and divisions... We should know this!!!
  } 
  else if (node.type=="path") #Never the case
  {
    state <- seqdecomp(path)[1]
    node.ccol <- Xtract("c.col", nPar, default="white")
    
    if (node.ccol=="state") {
      node.ccol <- if (path=="e") root.col else cpal[which(state==alphabet)]
    }
    symbols(X, Y, circles=node.size/2, inches=inches, add=TRUE, 
            fg=if (all(pruned)) pruned.col else par("col"), bg=node.ccol)
    ## State
    text(X, Y, state, xpd = TRUE, 
         ## srt = lab.srt, pos=lab.pos, offset=lab.offset,
         cex = lab.cex, col = lab.col, font = lab.font)
  }
  
  ## The node label: plotting the text in the square
  
  text(X, Y, nodeText, xpd = TRUE, srt = lab.srt, pos=lab.pos, offset=lab.offset,
       cex = lab.cex, col = lab.col, font = lab.font)
  
}

environment(plotNodeModified)=asNamespace("PST")

