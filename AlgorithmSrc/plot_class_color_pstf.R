####################################################################
# A function that plots the whole tree
# calls: plot_class_color_pstr
# Source: PST library
# Modified by: Gioia Wehbe
# last updated: 22-02-15
####################################################################


plot_class_color_pstf<-function (x, y, ...) 
{
  .local <- function (x, y = missing,title="", max.level = NULL, nodePar = list(), 
                      edgePar = list(), axis = FALSE, xlab = NA, ylab = if (axis) { "L"}
                      else {NA}, horiz = FALSE, xlim = NULL, ylim = NULL, withlegend = TRUE, 
                      ltext = NULL, cex.legend = 1, use.layout = withlegend != 
                        FALSE, legend.prop = NA, ...) 
  {
    if (nrow(x@cdata) > 0) 
    {
      ccol <- cpal(x@cdata)
      cnames <- alphabet(x@cdata)
      if (attr(x@cdata, "nr") %in% names(x[[2]])) 
      {
        ccol <- c(ccol, attr(x@cdata, "missing.color"))
        cnames <- c(cnames, attr(x@cdata, "nr"))
      }
      names(ccol) <- cnames
      if (!"stcol" %in% names(edgePar)) 
      {
        edgePar[["stcol"]] <- ccol
      }
      if (!"c.cpal" %in% names(nodePar)) {
        nodePar[["c.cpal"]] <- ccol
      } 
    }
    groups=levels(x@group)
    x <- as.pstree(x, max.level = max.level)
    plot_class_color_pstr(x, y = missing, title=title, max.level = max.level, groups=groups, nodePar = nodePar, 
         edgePar = edgePar, axis = axis, xlab = xlab, ylab = ylab, 
         horiz = horiz, xlim = xlim, ylim = ylim, withlegend = withlegend, 
         ltext = ltext, cex.legend = cex.legend, ...)
  }
  .local(x, y, ...)
}
environment(plot_class_color_pstf)=asNamespace("PST")

