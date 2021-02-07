#' pca.screeplot
#' 
#' Plots the percentage of explained variance against the number of the principal component.
#' 
#' @param data      The data used to build the screeplot
#' @param center    Center data before building screeplot
#' @param scale     Scale data before building screeplot
#' @param type      Type of screeplot (1 = dots and lines, 2 = barplot)
#' 
#' @examples
#'      data(iris)
#'      pca.screeplot(iris[,-5], scale = T)
#' 
#' @export
pca.screeplot <- function(data, center = T, scale = F, type = 1) {
  y <- scale(data, center = center, scale = scale)
  svd <- svd(y)
  
  if(type == 2) {
    graphics::barplot( svd$d^2/sum(svd$d^2)*100, names = c(1:ncol(data)), ylab="Percent variability explained", 
            xlab = "Principal Components", ylim = c(0,100))
  } else {
    graphics::plot( svd$d^2/sum(svd$d^2)*100, ylab="Percent variability explained", 
         xlab = "Principal Components", ylim = c(0,100), type = "b")
  }
}


