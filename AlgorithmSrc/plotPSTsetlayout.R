######################################################################################
# A function that sets the layout for TraMineR plots (including position of legend)
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 31-01-15
######################################################################################

PSTsetlayoutModified <- function(nplot, prows, pcols, withlegend, axes,legend.prop=NA) 
  {                             #  1 ,  NA   ,  NA  ,   TRUE    ,  all,     NA
  
  ## Backward compatibility
  if (withlegend==TRUE) withlegend <- "auto"
  
  if (is.na(pcols)) pcols <- min(nplot,2)   #set pcols = 1
  if (is.na(prows)) prows <- ceiling(nplot/pcols)   #set prows = 1
  
  ## Defining initial layout matrix
  pheight <- 1
  pwidth <- 1
  widths <- rep(pwidth/pcols,pcols)     #widths = c(1)
  heights <- rep(pheight/prows,prows)   #heights = c(1)
  
  layrow <- prows   #1
  laycol <- pcols   #1
  laymat <- matrix(1:(layrow*laycol), nrow=layrow, ncol=laycol, byrow=TRUE)   #a matrix of size 1x1 with value 1
  
  axisp <- 0
  
  legpos=NULL
  freecells <- (prows*pcols)-nplot    #0
  
  ## =========================
  ## Positioning of the legend
  ## =========================
  if (withlegend=="auto") 
    {
    if (freecells==0) #it is true
      {
      if (is.na(legend.prop))   #it is NA
#         legend.prop <- 0.15
        
        legend.prop <- 0.25
      layrow <- layrow+1  #layrow = 2
      
      pheight <- pheight-legend.prop    #1-0.15=0.85
      heights <- rep(pheight/prows,prows)   #heights=c(0.85)
      heights <- c(heights,legend.prop)   #heights=c(0.85,0.15)
      
      widths <- rep(pwidth/laycol,laycol)     #widths=c(1)
      
      legpos="bottom"
      
      ## Adding one row in the layout matrix for the legend
      laymat <- rbind(laymat, rep(nplot+1,ncol(laymat)))    # add the following row: c(2) so laymat is a matrix of 2 rows 1 column. first row value = 1, second row value = 2
    }                                 #2        #1
    else  #freecells is always = 0 so does not enter the else
      {
      legpos="center"
      heights <- rep(pheight/prows,prows)
      widths <- rep(pwidth/laycol,laycol)
    }
  }
  else if (withlegend=="right")   #whith legend is never = right unless I want to change the assignment if true
    {
    if (is.na(legend.prop)) legend.prop <- 0.25
    laycol <- laycol+1
    pwidth <- pwidth-legend.prop
    legpos="center"
    widths <- rep(pwidth/pcols,pcols)
    widths <- c(widths, legend.prop)
    heights <- rep(pheight/prows,prows)
    
    ## Adding one row in the layout matrix for the legend
    laymat <- cbind(laymat, rep(nplot+1,nrow(laymat)))
  }	
  
  ## if (axes %in% c("all","bottom")) axisp <- 1
  
  ## On which plots the axes will appear
  if (axes=="bottom") #axes = all so never enters if
    {
    for (nc in 1:ncol(laymat))
      axisp <- c(axisp, max(laymat[laymat[,nc]<=nplot,nc]))
  } 
  else if (axes=="all") axisp <- 1:nplot    #axisp = c(1)
  
  
  ## Returning a list with layout settings
  laylist <- list(laymat=laymat, widths=widths, heights=heights, axisp=axisp, legpos=legpos)
                      #1              C(1)        c(0.85,0.15)      c(1)          bottom
                      #2
  return(laylist)
} 


environment(PSTsetlayoutModified)=asNamespace("PST")


