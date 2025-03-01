#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------
plot_matrix_image = function(X, main = "Matrix",
                             labels = FALSE, 
                             save_plot = FALSE,
                             save_dir = ".",
                             file_name = "matrix_plot.png") {
  # If saving is requested, construct the full file path and open the device.
  if (save_plot) {
    file_path <- file.path(save_dir, file_name)
    png(filename = file_path, width = 800, height = 600)  # Adjust dimensions if needed.
  }
  
  # Plot the matrix.
  image(t(X)[, nrow(X):1], xaxt = "n", yaxt = "n", main = main)
  
  # Optionally add labels.
  if (labels) {
    axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
    axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
  }
  
  # If a device was opened, close it to save the plot.
  if (save_plot) {
    dev.off()
  }
}


# Funktion zum Plotten der Kovarianzmatrix
plot_matrix <- function(grid, sigma, phi) {
  # Erstellen der Kovarianzmatrix

  true_covariance <- cov_exponential(grid, sigma, phi, method = "difference")
  
  # Umwandeln der Matrix in ein langes Format für ggplot2
  cov_matrix_melted <- melt(true_covariance)
  
  # Plotten der Kovarianzmatrix
  ggplot(cov_matrix_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "x-Koordinate", y = "y-Koordinate", fill = "Wert", 
         title = paste("Wahre Kovarianzmatrix\nsigma =", sigma, ", phi1=phi2 =", phi)) +
    theme_minimal()
}

# Funktion zum Plotten der Kovarianzmatrix
plot_matrix_modi <- function(grid, sigma = 1, alpha1 = 1, alpha2 = 1,lambda1 =1, 
                             lambda2 = 1, beta = 0, test_sep = FALSE) {
  # Compute covariance matrix using ModifiedExponentialCovariance function
  cov_matrix <- ModifiedExponentialCovariance(grid, sigma, alpha1,
                                              alpha2, lambda1, lambda2, beta,
                                              test_sep)
  
  # Check if the result is a list (i.e., when test_sep = TRUE)
  if (is.list(cov_matrix)) {
    cov_matrix <- cov_matrix$covariance
  }
  
  # Melt the covariance matrix for ggplot2
  cov_melt <- melt(cov_matrix)
  
  # Create the plot
  ggplot(data = cov_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(midpoint = median(cov_melt$value),
                         low = "blue", mid = "white", high = "red",
                         space = "Lab", name="Covariance") +
    theme_minimal() +
    labs(x = "X", y = "Y", title = "Covariance Matrix") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
