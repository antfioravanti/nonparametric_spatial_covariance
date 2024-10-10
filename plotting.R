#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

plot_matrix = function(X, main = "Matrix", labels = F){
  # Plot Matrix by keeping fixed the cell positions
  
  image(t(X)[, nrow(X):1], xaxt = "n", yaxt = "n", main = main)
  
  if(labels == T){
    # Add x-axis labels
    axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
    
    # Add y-axis labels
    axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
  }
}
