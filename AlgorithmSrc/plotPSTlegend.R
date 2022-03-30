######################################################################################
# Plots the legend
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 31-11-14
######################################################################################


plotPSTlegend <- function(pos, text, colors, cex=1) 
  {
  nbstat <- length(text)    #num of classes = 3
  ## Computing some parameters for the legend's plotting
  if (pos=="bottom")  #always the case
    {     #Alignment of ledgend 
    if (nbstat > 6)   
      nbcol <- 6      #ledgend columns is the number of groups unless > 6 goups then fix to 6 and go to next row
    else
      nbcol <- nbstat
    leg.ncol <- nbcol
  }
  else
    leg.ncol <- 1   #Gioia
  ## leg.inset <- -0.2 + ((2-leg.ncol)*0.025)
  ## Setting graphical parameters while saving them in savepar
  savepar <- par(mar = c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
  ## Restoring graphical parameters
  on.exit(par(savepar))
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  ## legend(position, fill = cpal, legend = ltext, cex = fontsize)
  legend(pos,
         ## inset=c(0,leg.inset),         #Gioia
         legend=text,
         fill=colors,
         ncol=leg.ncol,
         bty="o",
         cex=cex)
}

environment(plotPSTlegend)=asNamespace("PST")


