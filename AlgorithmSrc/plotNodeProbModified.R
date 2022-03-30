#########################################################################
# Plots the nodes' square that represents the probabilities
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 31-01-15
#########################################################################


plotNodeProbModified<-function(x0, y0, x1, y1, prob,n=NULL, seglist, state, cpal, pruned, index, horiz=TRUE, 
         axes=c("no", "no"), bgcol="grey95", pruned.col="red", cex.axes=0.6, by.state=FALSE, type="b", 
         lwd=par("lwd"), x.prop=FALSE, frame=TRUE, ygrid=by.state, cm.pal="blue", xsrt=0, ysrt=0) 
{
  if (getOption("verbose")) 
  {
    cat("    [-] plotNodeProbModified: x0=", x0, ", y0=", y0, ", x1=", x1, ", y1=", y1, "\n")
  }

  A <- colnames(prob) #A is the alphabet or states
  
  if (!horiz) #Fixing the rectangle coordinates for the vertical display
  {
    x0.tmp <- x0
    x0 <- y0
    y0 <- x0.tmp
    x1.tmp <- x1
    x1 <- y1
    y1 <- x1.tmp
  }  
  
  xsize <- x1-x0    #Size of the horizontal sides of the rectangle
  ysize <- abs(y1-y0)     #Size of the vertical sides of the rectangle
  
  #seglist is the index slot of a node I think
  nbseg <- nrow(seglist)
  seg.lab <- rownames(seglist)
  has.segment <- rownames(seglist) %in% rownames(index)	
  rel_n=index[,"rel_n"]   #Gioia
  
  ## Filling the plotting area with the background color: simply draws a grey rectangle 
  rect(x0, y0, x1, y1, col=bgcol)
  
  pdim <- xsize/nbseg
  xleft <- x0
  
  if (type=="a") #not the case
  { 
    by.state <- TRUE 
    ygrid <- TRUE
  } 
  else if (type=="cm") #not the case
  {
    by.state <- TRUE 
    ygrid <- FALSE
    
    if (cm.pal=="blue") {
      ## cm.map <- rev(sequential_hcl(10, h = 210, c = c(80, 30), l = c(20, 80), power = 1.5))
      cm.map <- c("#9BCFDA", "#91C9D6", "#7FBFCD", "#65B2C1", "#40A3B4", "#0093A5", "#008195",
                  "#006F85", "#005D75", "#004E69")
    } else if (cm.pal=="heat") {
      cm.map <- c("#E2E6BD", "#E5DB80","#E4C96A","#E0B458","#D89E4B","#CE8642","#C16D3E",
                  "#B1523C", "#A0353C", "#8E063B")
    } else if (cm.pal=="heat20") { 
      cm.map <- c("#E2E6BD","#E4E38E","#E5DC81","#E5D476","#E4CB6C","#E3C263","#E1B85B",
                  "#DEAD53","#DAA24D","#D69748","#D18C44","#CB8141","#C5753F","#BE693E","#B75C3D",
                  "#B04F3C","#A8423C","#9F333C","#97223C","#8E063B")
    }
  }
  
  if (by.state) #by.state is always default = FALSE
  {
    ## stsep <- (0.10*ysize)/length(A)
    stsep <- 0
    stsize <- (ysize-((length(A)-1)*stsep))/length(A)
    ytmp <- y0
  } 
  else 
  { 
    stsize <- ysize   #stsize=ysize=size of the vertical side of the rectangle
  }
  
  if (ygrid) #ygrid is always default = FALSE
  {	
    for (s in 1:(length(A)-1)) 
    {
      ytmp <- ytmp-stsize
      segments(x0, ytmp, x1, ytmp, col="grey50")
    }
  }
  
  ## NOW HERE IS THE BULCK OF THE WHOLE THING .... HERE WE SHOULD MODIFY... we should divide the rectangle only based on the group and not the state!!!
  for (g in 1:nbseg) #for every group
  {
    xright <- xleft+pdim
    ytmp <- y0
    ##########################################Gioia
    if (has.segment[g])   #if the current node belongs to this group...
    {
      idseg <- rownames(seglist)[g]
      rel_n=index[idseg,"rel_n"]
      ybot=ytmp-(stsize*rel_n)
      if (type=="b") #Always the case!!
      { 
        rect(xleft, ytmp, xright, ybot, col=cpal[g], border=NA)
        ytmp <- if (by.state) { ytmp-stsize } else { ybot }   #ytemp is always == ybot
      }
      
    }
    xleft <- xright 
  }	#end of for (g in 1:nbseg)
  
  
  ## Frame surrounding the area
  if (frame) { rect(x0, y0, x1, y1) }
  
  ## Plotting the axes
  ## x axis
  if (axes[1]=="bottom") {
    axe.offset <- if (nbseg>1 && any(pruned, na.rm=TRUE)) { ysize*0.3 } else {0}
    segments(x0+(pdim/2), y0, x0+(pdim/2), y0+(0.10*ysize))
    segments(x1-(pdim/2), y0, x1-(pdim/2), y0+(0.10*ysize))
    text(x=c(x0+(pdim/2),x1-(pdim/2)), y=y0+(0.25*ysize), labels=c(1,nbseg), cex=cex.axes, srt=xsrt)
  } else if (axes[1]=="top") {
    segments(x0+(pdim/2), y1, x0+(pdim/2), y1-(0.10*ysize))
    segments(x1-(pdim/2), y1, x1-(pdim/2), y1-(0.10*ysize))
    text(x=c(x0+(pdim/2),x1-(pdim/2)), y=y1-(0.25*ysize), labels=c(1, nbseg), cex=cex.axes, srt=xsrt)
  }
  
  ## y axis
  if (!by.state) {
    if (axes[2]=="left") {
      ## segments(x0-(0.1*xsize), y0, x0-(0.1*xsize), y1)
      segments(x0, y0, x0-(0.10*xsize), y0)
      segments(x0, y1, x0-(0.10*xsize), y1)
      segments(x0, (y0+y1)/2, x0-(0.05*xsize), (y0+y1)/2)
      text(x=c(x0-(0.25*xsize),x0-(0.25*xsize)), y=c(y0, y1), labels=c(0,1), cex=cex.axes, srt=ysrt)
    } else if (axes[2]=="right") {
      segments(x1+(0.1*xsize), y0, x1+(0.1*xsize), y1)
      segments(x1+(0.1*xsize), y0, x1+(0.15*xsize), y0)
      segments(x1+(0.1*xsize), y1, x1+(0.15*xsize), y1)
      text(x=c(x1+(0.25*xsize),x1+(0.25*xsize)), y=c(y0, y1), labels=c(0,1), cex=cex.axes, srt=ysrt)
    }
  }
  
  ## A bar showing the pruned and unpruned nodes
  if (nbseg>1 && any(pruned, na.rm=TRUE)) {
    rect(x0, y0+(ysize*0.1), x1, y0+(ysize*0.3), col=bgcol, border=NA)
    xleft <- x0
    for (g in 1:nbseg) {
      xright <- xleft+pdim
      if (has.segment[g] %in% rownames(prob)) {
        if (pruned[has.segment[g]]) {
          rect(xleft, y0+(ysize*0.1), xright, y0+(ysize*0.3),
               col=pruned.col, border=NA)
        } else if (!pruned[has.segment[g]]) {
          rect(xleft, y0+(ysize*0.1), xright, y0+(ysize*0.3), 
               col="green", border=NA)
        }
      }
      xleft <- xright
    }
    rect(x0, y0+(ysize*0.1), x1, y0+(ysize*0.3))
  } else if (nbseg==1 && pruned) {
    segments(x0, y0, x1, y1, col = pruned.col)
  }
}




environment(plotNodeProbModified)=asNamespace("PST")






