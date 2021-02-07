#' pca.biplot
#' 
#' A customizable correlational biplot function. PC's are scaled to unit variance.
#' 
#' @param data            The data used to build the biplot
#' @param center          Column centering the data (default = \code{TRUE}). If center is \code{TRUE} then centering is done by subtracting the column means (omitting \code{NA}'s) of x from their corresponding columns, and if center is \code{FALSE}, no centering is done.
#' @param scale           Scaling the data (default = \code{FALSE}). If scale is \code{TRUE} then scaling is done by dividing the (centered) columns of x by their standard deviations if center is \code{TRUE}, and the root mean square otherwise. If scale is \code{FALSE}, no scaling is done.
#' @param pcx             Which principal component to plot on the X-axis (default = 1)
#' @param pcy             Which principal component to plot on the Y-axis (default = 2)
#' @param groups          A vector from which the unique values are used to color the observations
#' @param pch             Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link[graphics]{points}} for possible values and their interpretation.
#' @param col             Vector of colors, to color observations. Needs to be equal to the number of unique values of \code{groups} variable. If no \code{groups} is defined then all points will have the defined color.
#' @param ellipse         Whether to add confidence ellipses to the different groups
#' @param ellipse.alpha   Alpha transparency value for the fill color (0 = transparent, 1 = opaque)
#' @param ellipse.conf    Confidence limit for ellipses, e.g. 0.68. If given, the corresponding sd or se is multiplied with the corresponding value found from the Chi-squared distribution with 2df.
#' @param ellipse.lwd     Line width of ellipse border
#' @param arrow.options   List to define options for the variable arrows (for all options see \code{\link[graphics]{arrows}} ). Needs to be defined as a list, e.g. \code{arrows.options = list(angle = 0.15, col = "blue")}, errors could occur if a vector is used.
#' @param label.options   List to define options for the variable labels (for all options see \code{\link[graphics]{text}} ). Needs to be defined as a list, e.g. \code{label.options = list(cex = 1.2, col = "red")}, errors could occur if a vector is used.
#' @param abbrev          Whether to abbreviate the variable names
#' @param circle          Whether to add a unit circle. Only applies if data is scaled!
#' @param circle.options  List to define options for the unit circle (for all options see \code{\link[graphics]{polygon}} ). Needs to be defined as a list, e.g. \code{circle.options = list(lty = "dashed", col = "grey")}, errors could occur if a vector is used.
#' @param legend          Whether to add a legend
#' @param legend.pos      Position of the legend
#' @param ...             Other graphical parameters, see \code{\link[graphics]{plot}} for more details.
#' 
#' @examples 
#'   data(iris)
#'   pca.biplot(iris[1:4])
#'   pca.biplot(iris[,-5], col = c("#00AFBB", "#E7B800", "#FC4E07"), groups = iris$Species, 
#'        legend = F, arrow.options = list(col="blue", angle = 20, length = 0.15))
#'   pca.biplot(iris[,-5], scale = T, col = c("#00AFBB", "#E7B800", "#FC4E07"), groups = iris$Species, 
#'        ellipse = T, pch = 16, legend.pos = 'topright')
#' @export
pca.biplot <- function(data, center = TRUE, scale = FALSE, pcx = 1, pcy = 2, groups = 1,
                       col = NULL, pch = 1,
                       ellipse = FALSE, ellipse.alpha = 0.2, ellipse.conf = 0.68, ellipse.lwd = 2,
                       arrow.options,
                       label.options, abbrev = FALSE,
                       circle = TRUE, circle.options,
                       legend = TRUE, legend.pos = "bottomleft",
                       ...
){
  data <- as.data.frame( lapply(data, as.numeric) )     # Make sure everything is numeric
  y <- scale(data, center = center, scale = scale)    # center and scale
  svd <- svd(y)                                     # SVD
  diag.d <- diag(svd$d)
  
  if(is.null(col)){
    ellipse.col <- grDevices::palette()
    col <- as.factor(groups)
  } else if(length(col) != length(unique(groups))) {
    # Check if length col vector is equal to number of unique groups
    stop("col not defined correctly, number of colors should be equal to the number of unique groups")
  } else {
    ellipse.col <- col
    col <- col[as.factor(groups)]
  }
  
  # standardized PC scores
  sqrtn <- sqrt(nrow(svd$u) - 1)
  x.scores <- (svd$u * sqrtn)[,pcx]
  y.scores <- (svd$u * sqrtn)[,pcy]
  size <- max( abs( c( range(x.scores), range(y.scores) ) ) )
  
  # Plot PC's standardized to unit variance
  plot.new()
  graphics::plot(x = x.scores, y = y.scores, xlim = c(-size, size), ylim = c(-size, size),
                 xlab = paste("PC ", pcx, " (", round((svd$d^2/sum(svd$d^2)*100)[pcx], digits = 1), "% explained var.)", sep = "" ),
                 ylab = paste("PC ", pcy, " (", round((svd$d^2/sum(svd$d^2)*100)[pcy], digits = 1), "% explained var.)", sep = "" ),
                 col = col, pch = pch, ...)
  
  # Add confidence ellipse to PCA plot, only if levels are defined using colors
  if(ellipse && length(unique(col)) > 1 && exists("groups") && length(unique(col)) == length(unique(groups)) ) {
    for(i in 1:nlevels(as.factor(groups))) {
      vegan::ordiellipse((svd$u * sqrtn), groups, draw = "polygon", conf = ellipse.conf, lwd = ellipse.lwd, col = ellipse.col[i], show.groups = levels(as.factor(groups))[i], border = ellipse.col[i], alpha = ellipse.alpha)
    }
  }
  
  # PC loadings
  loadings <- (svd$v%*%diag.d)/sqrtn
  x.loadings <- loadings[,pcx]
  y.loadings <- loadings[,pcy]
  loadings.size <- max( abs( c( range(x.loadings), range(y.loadings) ) ) )
  
  # Add new plot window to plot PC loadings on. 
  # Plotting them on the same window as observations makes loadings quite small on the plot
  graphics::par(new = TRUE, las = 1)
  if(scale == F) {
    graphics::plot.window(xlim = c(-loadings.size, loadings.size), ylim = c(-loadings.size, loadings.size))
  } else if (scale == T) {
    graphics::plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
  }
  
  # Add variable arrows to new plot window
  if(missing(arrow.options)) {
    graphics::arrows(x0 = 0, x1 = loadings[,pcx], y0 = 0, y1 = loadings[,pcy], col = "black", length = 0.08, angle = 30, lwd = 1)
  } else {
    do.call(graphics::arrows, c(list(x0 = 0, x1 = loadings[,pcx], y0 = 0, y1 = loadings[,pcy]), arrow.options))
  }
  
  # Abbreviate Variable Names ?
  if(abbrev) {
    label.name <- as.character( abbreviate( dimnames(data)[[2]] ) )
  } else {
    label.name <- as.character( dimnames(data)[[2]] )
  }
  
  # Add variable labels to arrows
  if(missing(label.options)) {
    graphics::text(loadings[,pcx]-.05, loadings[,pcy]-.05, labels = label.name , cex = 0.8)
  } else {
    do.call(graphics::text, c(list(loadings[,pcx]-.05, loadings[,pcy]+.05, labels = label.name), label.options))
  }
  
  # Add unit circle (only if data is scaled)
  if(circle && scale == TRUE){
    ucircle = cbind(cos((0:360)/180*pi), sin((0:360)/180*pi))
    if(missing(circle.options)) {
      graphics::polygon(ucircle)
    } else {
      do.call(graphics::polygon, c(list(ucircle), circle.options))
    }
  }
  
  # Add a legend explaining the colors
  if(legend == TRUE && length(groups) > 1) {
    graphics::legend(legend.pos, legend = unique(groups), col = unique(col), pch = pch)
  }
  
}
