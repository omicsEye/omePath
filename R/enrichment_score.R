# simulate two samples
enrichment_score <- function(){
  a <- stats::rnorm(100)
  b <- stats::rnorm(100, 2)
  
  # define limits of a common grid, adding a buffer so that tails aren't cut off
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  da <- stats::density(a, from=lower, to=upper)
  db <- stats::density(b, from=lower, to=upper)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  total <- sfsmisc::integrate.xy(d$x, d$a) + sfsmisc::integrate.xy(d$x, d$b)
  intersection <- sfsmisc::integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  overlap <- 2 * intersection / total
  
  plot(d$x, d$a)
  plot(d$x, d$b)
  
  ### Approach 2
  # https://stackoverflow.com/questions/41914257/calculate-area-of-overlapping-density-plot-by-ggplot-using-r
  
  ##  Create the two density functions and display
  FDensity = stats::approxfun(stats::density(df$weight[df$sex=="F"], from=40, to=80))
  MDensity = stats::approxfun(stats::density(df$weight[df$sex=="M"], from=40, to=80))
  plot(FDensity, xlim=c(40,80), ylab="Density")
  graphics::curve(MDensity, add=TRUE)
  
  ## Solve for the intersection and plot to confirm
  FminusM = function(x) { FDensity(x) - MDensity(x) }
  Intersect = stats::uniroot(FminusM, c(40, 80))$root
  graphics::points(Intersect, FDensity(Intersect), pch=20, col="red")
  
  
  #Now we can just integrate to get the area of the overlap.
  
  stats::integrate(MDensity, 40,Intersect)$value + 
    stats::integrate(FDensity, Intersect, 80)$value
}
