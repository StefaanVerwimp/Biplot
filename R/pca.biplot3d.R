#' pca.biplot3d
#' 
#' An interactive 3D biplot. PC's are scaled to unit variance.
#' 
#' @param data            The data used to build a biplot
#' @param center          Whether the variables should be shifted to be zero centered (default = TRUE)
#' @param scale           Whether the variables should be scaled to have unit variance (default = FALSE)
#' @param pcx             Which principal component to plot on the X-axis (default = 1)
#' @param pcy             Which principal component to plot on the Y-axis (default = 2)
#' @param pcz             Which principal component to plot on the Z-axis (default = 3)
#' @param col             Change color of observations, only if no group is defined
#' @param groups          Color observations based on the unique values of this group
#' @param bg              Background color
#' @param type            Type of plot. Options: "sphere", "point", "text" (default = "sphere")
#' @param r               Vector or single value defining the sphere radius/radii
#' @param aspect3d        Vector used to set apparent ratios of the x, y, and z axes of the current bounding box
#' @param abbrev          Whether to abbreviate the variable names
#' @param arrow.col       Color of the variable arrow
#' @param label.col       Color of the variable label
#' @param label.cex       Size of the variable label
#' @param simple.axis     Use axes3d instead of decorate3d if TRUE
#' @param box             Whether to draw a box (only if using decorate3d)
#' @param grid            Adds a reference grid
#' @param ellipse         Whether to add confidence ellipses to groups (to display this on the plot, ellipse.shade and/or ellipse.wire must be TRUE)
#' @param ellipse.conf    Confidence limit for ellipses
#' @param ellipse.alpha   Alpha transparency value for the ellipse (0 = transparent, 1 = opaque, default = 0.2)
#' @param ellipse.shade   Draws ellipse as a surface
#' @param ellipse.wire    Draws ellipse as line segments
#' 
#' @param ...   Material properties (see rgl.material).
#' 
#' @examples 
#'      data(iris)
#'      pca.biplot3d(iris[,-5], scale = T)
#'      pca.biplot3d(iris[,-5], groups = iris$Species, scale = T, box = T, grid = T, ellipse = T)
#' 
#' @export
pca.biplot3d <- function(data, center = TRUE, scale = FALSE, pcx = 1, pcy = 2, pcz = 3, col = "black", groups = NULL,
                      bg = "white", type = "sphere", r = 0.05, aspect3d = c(1, 1, 1), abbrev = FALSE,
                      arrow.col = "black", label.col = "darkblue", label.cex = 1,
                      simple.axis = FALSE, box = FALSE, grid = FALSE,
                      ellipse = FALSE, ellipse.conf = 0.68, ellipse.alpha = 0.2, ellipse.shade = TRUE, ellipse.wire = FALSE,
                      ...
){
  data <- as.data.frame( lapply(data, as.numeric) )     # Make sure everything is numeric
  y <- scale(data, center = center, scale = scale)    # center and scale
  svd <- svd(y)                                     # SVD
  
  # PC scores
  sqrtn <- sqrt(nrow(svd$u) - 1)
  x.scores <- (svd$u * sqrtn)[,pcx]
  y.scores <- (svd$u * sqrtn)[,pcy]
  z.scores <- (svd$u * sqrtn)[,pcz]
  size <- max( abs( c( range(x.scores), range(y.scores) , range(z.scores)) ) )
  
  # Colors
  if(!is.null(groups)) {
    col <- get_colors(groups)
  }
  
  rgl::rgl.open()
  rgl::rgl.bg(color = bg)
  
  if (type == "points") {
    rgl::points3d(x.scores, y.scores, z.scores,  color = col, ...)
  } else if (type == "text") {
    rgl::text3d(x.scores, y.scores, z.scores, texts=as.character(c(1:length(x.scores))), col = col, ...)
  } else {
    rgl::rgl.spheres(x.scores, y.scores, z.scores, r = r, color = col, aspect = F, ...)
  }
  
  # Adds ellipse to observations
  if(ellipse && !is.null(groups)) {
    lvls <- levels(groups)
    group.col <- unique(col)
    for(i in 1:length(lvls)) {
      group <- lvls[i]
      selected <- groups == group
      xx <- x.scores[selected]; yy <- y.scores[selected]; zz <- z.scores[selected]
      ellips <- rgl::ellipse3d(cov(cbind(xx, yy, zz)), centre = c(mean(xx), mean(yy), mean(zz)), level = ellipse.conf)
      if(ellipse.shade) rgl::shade3d(ellips, col = group.col[i], alpha = ellipse.alpha, lit = FALSE)
      if(ellipse.wire) rgl::wire3d(ellips, col = group.col[i], alpha = ellipse.alpha,  lit = FALSE)
      rgl::texts3d(mean(xx),mean(yy), mean(zz), text = group, col= group.col[i], cex = label.cex)
    }
  }
  
  # Sets the apparent ratios of the x, y, and z axes of the current bounding box.
  rgl::aspect3d(aspect3d)
  
  # PC loadings
  diag.d <- diag(svd$d)
  loadings <- (svd$v%*%diag.d)/sqrtn
  
  # Plot variable arrows
  for ( i in 1:nrow(loadings) ) {
    rgl::segments3d(suppressWarnings(rbind(matrix(0, nc = 3), loadings[i,])), col = arrow.col)
  }
  
  # Abbreviate variable labels?
  if(abbrev) {
    label.name <- as.character( abbreviate( dimnames(data)[[2]] ) )
  } else {
    label.name <- as.character( dimnames(data)[[2]] )
  }
  # Variable labels
  rgl::text3d(loadings-0.05, texts = label.name, col = label.col, cex = label.cex)
  
  # Add axis
  if(simple.axis) {
    rgl::axes3d(c('x.scores', 'y.scores', 'z.scores'), col = "black") 
  } else {
    rgl::decorate3d(xlab = "", ylab = "", zlab = "", box = box)
  }
  
  # Label axis
  rgl::title3d(xlab = paste("PC ", pcx, " (", round((svd$d^2/sum(svd$d^2)*100)[pcx], digits = 1), "%)", sep = "" ),
               ylab = paste("PC ", pcy, " (", round((svd$d^2/sum(svd$d^2)*100)[pcy], digits = 1), "%)", sep = "" ),
               zlab = paste("PC ", pcz, " (", round((svd$d^2/sum(svd$d^2)*100)[pcz], digits = 1), "%)", sep = "" ),
               color = "black")
  
  # Add grid?
  if(grid) rgl::grid3d(c("x", "y", "z"))
  
}

get_colors <- function(groups, group.col = grDevices::palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

