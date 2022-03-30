#########################################################################
# A function that plots the pstr node (subtree)
# calls: PSTsetlayoutModified, plotTreeModified,plotPSTlegend
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 22-02-15
#########################################################################


plot_class_color_pstr<-function (x, y, ...) 
{
  .local <- function (x, y = missing, title="", max.level = NULL, groups=NULL, nodePar = list(), 
                      edgePar = list(), axis = FALSE, xlab = NA, ylab = if (axis) {"L"}
                      else {NA}, horiz = FALSE, xlim = NULL, ylim = NULL, withlegend = TRUE, 
                      ltext = NULL, cex.legend = 1, use.layout = withlegend != 
                        FALSE, legend.prop = NA, ...) 
  {
    #Gioia: represents color of each class instead of each state. updated to represent the color of each class from 
    #the pstree function
    cpal <- x@cpal    
    Lmar <- if (axis) {4}
    else {2}
    if (horiz) {
      par(mar = c(Lmar, 2, 4, 2))
    }
    else {
      par(mar = c(2, Lmar, 4, 2))
    }
    oolist <- list(...)     
    cex <- if ("cex" %in% names(oolist)) {    #cex has to do with the size of the text 
      oolist[["cex"]]
    }
    else {  #cex always = 1
      1
    }
    xaxt <- if ("xaxt" %in% names(oolist)) 
      {    #what is xaxt DOESNOT ENTER THIS IF
      oolist[["xaxt"]]
    }
    else {  #xaxt always = "n"
      "n"
    }
    yaxt <- if ("yaxt" %in% names(oolist)) {    #what is yaxt DOESNOT ENTER THIS IF
      oolist[["yaxt"]]
    }
    else {      #yaxt always = "n"
      "n"
    }
    if (use.layout) #use.layout is true if legend is true and legend is true by default
      {
      savepar <- par(no.readonly = TRUE)
      lout <- PSTsetlayoutModified(nplot = 1, prows = NA, pcols = NA, 
                            withlegend, axes = "all", legend.prop)  #this specifies the layout of the plot including legend position
      layout(lout$laymat, heights = lout$heights, widths = lout$widths)     #layout divides the device up into as many rows and columns as there are in matrix laymat, with the column-widths and the row-heights specified in the respective arguments.
                  #1      row heights: c(0.75,0.25)  col widths: c(1)          
                  #2      
      
      legpos <- lout$legpos     #bottom
    }
    else    #never the case
      {
      legpos <- NULL
    }
    
    
    seglist <- x@index         #index includes group of current node and position (not implimented yet)
    k <- x@order               #depth of the node / length of the substring it represents 
    stats <- summary(x, max.level = max.level, segmented = FALSE)
    if (missing(max.level) | is.null(max.level)) {
      max.level <- stats@depth        #maximum depth=4
      max.level=max.level+1
    }
    if (!"cpal" %in% names(nodePar)) #cpal is in node par 
      {
      nodePar[["cpal"]] <- attr(x, "cpal")
    }
    
    hgt <- max.level      #height of the tree
    mem.x <- stats@leaves    #number of leaves
    pin <- par("pin")       #The current plot dimensions, (width, height), in inches.
    node.type <- Xtract("node.type", nodePar, default = "prob")   #node.type is "prob"
    
    node.size <- Xtract("node.size", nodePar, default = min(0.6,((mem.x - 1)/mem.x) - 0.1))
#     node.size <- Xtract("node.size", nodePar, default = 0.6)
    if(mem.x==1)
      node.size <- Xtract("node.size", nodePar, default = 0.1)
    if (!"node.size" %in% names(nodePar)) {
      nodePar[["node.size"]] <- node.size
    }
    #gratio is the ratio btwn horizontal and vertical dimentions of a node
   
    gratio <- Xtract("gratio", nodePar, default = min((((hgt - k) + 1)/mem.x), 1))

    #what is leave.lh and leave.lw, leave height and leave width?
    leave.lh <- Xtract("leave.lh", edgePar, default = 0.1)    
    leave.lw <- Xtract("leave.lw", edgePar, default = node.size)

    
    yTop <- k #depth of node
    x1 <- 1
    x2 <- mem.x       #number of leaves
    if (horiz) 
      {  #if horizontal
      xl. <- c(x1 - ((node.size/2) * gratio), x2 + ((node.size/2) *gratio))     #x-limits
      yl. <- c(k - (node.size/2), hgt + (node.size/2) +  leave.lh)              #y-limits
    }    #Not the case
    else 
      {#if vertical
      ym <- if (node.type == "prob") 
        {
        0.5
      }
      else 
      {
        (node.size/2) * gratio
      }
      xl. <- c(x1 - (node.size/2), x2 + (node.size/2))      #x-limits
      yl. <- c(k - ym, hgt + ((node.size/2) * gratio) +leave.lh)    #y-limits
    }
    yl. <- rev(yl.)   #reverse y-limit
    if (horiz) #if horizontal, switch the x and y of x&y limit, x&y axt and x&y label
    {
      tmp <- xl.
      xl. <- yl.    #xlimit = ylimit
      yl. <- rev(tmp) #ylimit = reverse of original xlimit
      tmp <- xaxt   #what is xaxt
      xaxt <- yaxt    #set xaxt to yaxt
      yaxt <- tmp   #set yaxt to xaxt
      tmp <- xlab
      xlab <- ylab  #similarly, switch xlab and ylab
      ylab <- tmp
    }
    if (missing(xlim) || is.null(xlim)) { #set xlim
      xlim <- xl.
    }
    if (missing(ylim) || is.null(ylim)) {    #set ylim
      ylim <- yl.
    }
    plot(0, main=paste(title,node.size,mem.x),xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
         frame.plot = FALSE, ylab = ylab, xaxt = xaxt, yaxt = yaxt, #it only fixes and sets up the window for plotting. 
         ...)
    if (horiz) 
    {
      nc <- par("pin")[2]/par("usr")[3]   #specify some kind of coordinate = 7.025463
    }
    else {
      nc <- NULL
    }
    if (axis) 
    {
      if (horiz) 
      {
        axis(1, at = 0:hgt)#set the axis
      }
      else {
        axis(2, at = 0:hgt)#set the axis
      }
    }
    
    
    plotTreeModified(x1, x2, x, seglist, nPar = nodePar, ePar = edgePar,    #gioia: this is the function to modify
                     horiz = horiz, gratio = gratio, max.level = max.level, 
                     cex = cex, nc = nc, cpal = cpal)       
    if (!is.null(legpos)) ###Legend position is not null (default)
    {
      if(length(groups)==0)
      {
        index=x@index
        groups=index[,'group']
      }
      if (is.null(ltext)) 
        ltext <- groups    #text of the legend (labels of the legend modified to be classes)
      plotPSTlegend(legpos, ltext, cpal, cex = cex.legend)
    }
    if (use.layout) {
      par(savepar)
    }
  }
  .local(x, y, ...)
}
environment(plot_class_color_pstr)=asNamespace("PST")
