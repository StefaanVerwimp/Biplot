#' pca.screeplot
#' 
#' Plots the percentage of explained variance against the number of the principal component.
#' 
#' @param data      The data used to build the screeplot
#' @param center    Column centering the data (default = `TRUE`). If center is `TRUE` then centering is done by subtracting the column means (omitting `NAs`) of x from their corresponding columns, and if center is `FALSE`, no centering is done.
#' @param scale     Scaling the data (default = `FALSE`). If scale is `TRUE` then scaling is done by dividing the (centered) columns of x by their standard deviations if center is `TRUE`, and the root mean square otherwise. If scale is `FALSE`, no scaling is done.
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


