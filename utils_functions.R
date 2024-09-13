# SCRIPT OF FUNCTIONS
#-------------------------------------------------------------------------------

if (!require(MASS)) install.packages("MASS"); library(MASS)

#-------------------------------------------------------------------------------

is_positive_semi_definite = function(matrix) {
  # Function to check that a matrix is positive semi definite
  eigenvalues <- eigen(matrix)$values
  return(all(eigenvalues >= 0))
}

#-------------------------------------------------------------------------------
# Analytical Covariances
#-------------------------------------------------------------------------------

cov_exponential = function(grid, sigma, phi, method = "euclidean") {
  # COMPUTES THE ANALYTICAL COVARIANCE MATRIX WITH EXPONENTIAL
  # Can either be based on eucleadian distance or on the absolute difference
  if(!sum(method == c("difference", "euclidean"))){
    print("Specify either euclidean or difference method")
  }

  if(method == "euclidean"){
    # Euclidean distance
    dist_matrix = as.matrix(dist(grid, method = method))
    
    covariance = sigma^2 * exp(-dist_matrix / phi)
    
  }else if (method == "difference"){
    # Absolute difference
    dist_matrix_x = abs(outer(grid$x,grid$x,"-"))
    dist_matrix_y = abs(outer(grid$y,grid$y,"-"))
    
    covariance = sigma^2*exp(-dist_matrix_x / phi)*exp(-dist_matrix_y / phi)
  }
  return(covariance)
}

# MODIFIED EXPONENTIAL COVARIANCE


# MODIFIED EXPONENTIAL
ModifiedExponentialCovariance = function(grid,
                                         sigma = 1,
                                         alpha1 = 1,
                                         alpha2 = 1,
                                         lambda1 = 1,
                                         lambda2 = 1,
                                         beta = 0,
                                         test_sep = F) {
  # Args
  #   alpha1 and alpha 2 control the smoothness in each dimension
  #   lambda1 and lamda 2 control the range in each dimension
  #   beta controls the separability
  #     beta = 0  separable
  #     beta = 1 non-separable
  
  
  # Check that grid has the correct column names "x" and "y"
  expected_names <- c("x", "y")
  
  if (!all(colnames(grid) == expected_names)) {
    stop("Error: The grid must have column names 'x' and 'y'.")
  }
  
  # Compute pairwise differences for each dimension
  h1_diff = outer(grid$x, grid$x, "-")  # Difference in first dimension
  h2_diff = outer(grid$y, grid$y, "-")  # Difference in second dimension
  
  cov = sigma^2 * exp(-(abs(h1_diff)^alpha1 / lambda1 +
                          abs(h2_diff)^alpha2 / lambda2 + 
                          beta * abs(h1_diff - h2_diff)))
  if(test_sep == T){
    cov_sep = sigma^2 * exp(-(abs(h1_diff)^alpha1 / lambda1)-beta*abs(h1_diff)) * 
      exp(-(abs(h2_diff)^alpha2 / lambda2) -beta*abs(h2_diff))
    
    norm_diff = norm(cov_sep - cov, type = "2")
    return(list(covariance = cov, norm_diff = norm_diff))
  }else{
    return(cov)
  }
}
#-------------------------------------------------------------------------------

# Funktion zum Plotten der Kovarianzmatrix
plot_matrix <- function(grid, sigma, phi) {
  # Erstellen der Kovarianzmatrix
  true_covariance <- cov_exponential(grid, sigma, phi, method = "euclidean")
  
  # Umwandeln der Matrix in ein langes Format für ggplot2
  cov_matrix_melted <- melt(true_covariance)
  
  # Plotten der Kovarianzmatrix
  ggplot(cov_matrix_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "Index 1", y = "Index 2", fill = "Covariance", 
         title = paste("Covariance Matrix\nsigma =", sigma, ", phi =", phi)) +
    theme_minimal()
}

plot_matrix_modi <- function(grid, sigma = 1, alpha1 = 1, alpha2 = 1, lambda1 = 1, 
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

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#  Kernels
#-------------------------------------------------------------------------------

gaussian_kernel = function(x, bandwidth) {
  # GAUSSIAN KERNEL FUNCTION
  return((1 / (sqrt(2 * pi) * bandwidth)) * exp(-0.5 * (x / bandwidth)^2))
}

epanechnikov_kernel = function(x, bandwidth) {
  # EPANECHNIKOV KERNEL FUNCTION
  k = 0.75 * (1 - (x / bandwidth)^2)
  k[abs(x) > bandwidth] = 0  # Set the kernel to 0 outside the bandwidth
  return(k)
}

rechteck_kernel <- function(x, bandwidth) {
  # RECHTECK KERNEL FUNKTION
  return(ifelse(abs(x / bandwidth) <= 0.5, 1 / bandwidth, 0))
}

triangle_kernel <- function(x, bandwidth) {
  # TRIANGLE KERNEL FUNKTION
  return(ifelse(abs(x / bandwidth) <= 1, (1 - abs(x / bandwidth)), 0))
}

cubic_b_spline_kernel <- function(x, bandwidth) {
  # Normalize the distance
  u <- x / bandwidth
  
  # Cubic B-spline kernel formula
  kernel_value <- 0.5 * exp(-abs(u) / sqrt(2)) * (sin(abs(u) / sqrt(2)) + pi / 4)
  
  return(kernel_value)
}

#-------------------------------------------------------------------------------
# Estimators 
#-------------------------------------------------------------------------------

weight_function <- function(t1,t2){
  log <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  ifelse(abs(log) <= 10, 1, 0)
}

# Vektorisierte Kernel-Kovarianzfunktion

kernel_cov = function(t1, t2, X, grid, bandwidth,
                      kernel_function="gaussian_kernel"){
  # COMPUTES THE COVARIANCE ESTIMATOR
  
  # Ensure the kernel_function is correctly specified
  kernel_fun = switch(
    kernel_function,
    gaussian_kernel = gaussian_kernel,
    epanechnikov_kernel = epanechnikov_kernel,
    rechteck_kernel = rechteck_kernel,
    triangle_kernel = triangle_kernel,
    cubic_b_spline_kernel = cubic_b_spline_kernel,
    stop("Specify either gaussian_kernel or epanechnikov_kernel
         or rechteck_kernel or triangle_kernel kernel 
         or or cubic_b_spline_kernel")
  )
  
  Xbar = mean(X) # Sample mean of the simulated values X
  demeaned_X = X - Xbar # X- Xbar
  var = mean((X - Xbar)^2)  # Varianz der Stichprobe
  
  # Lags between the two points
  lag_x = t1[1] - t2[1]
  lag_y = t1[2] - t2[2]
  
  # Berechnet die Differenzen zwischen allen x-Koordinaten der Punkte im Raster
  diff_x = outer(grid$x, grid$x, "-")
  # Berechnet die Differenzen zwischen allen y-Koordinaten der Punkte im Raster
  diff_y = outer(grid$y, grid$y, "-")
  
  # Apply the chosen kernel
  k_x = kernel_fun(lag_x - diff_x, bandwidth)
  k_y = kernel_fun(lag_y - diff_y, bandwidth)
  K_vals = k_x * k_y

  X_ij = outer(demeaned_X, demeaned_X, "*") 
  
  cov_vals = X_ij / var
  
  
  # Numerator
  numerator = sum(K_vals * cov_vals)
  # Denominator
  denominator = sum(K_vals)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  return(list(covariance = numerator / denominator, weights = K_vals))
  
}

#-------------------------------------------------------------------------------
# Vektorisierte Kernel-Kovarianzfunktion mit euklidischer Distanz

kernel_cov_euclidean <- function(t1, t2, X, grid, bandwidth,
                                 kernel_function="gaussian_kernel") {
  # COMPUTES THE COVARIANCE ESTIMATOR BASED ON EUCLIDEAN DISTANCE
  
  # Ensure the kernel_function is correctly specified
  kernel_fun = switch(
    kernel_function,
    gaussian_kernel = gaussian_kernel,
    epanechnikov_kernel = epanechnikov_kernel,
    rechteck_kernel = rechteck_kernel,
    triangle_kernel = triangle_kernel,
    cubic_b_spline_kernel = cubic_b_spline_kernel,
    stop("Specify either gaussian_kernel or epanechnikov_kernel
         or rechteck_kernel or triangle_kernel kernel 
         or or cubic_b_spline_kernel")
  )
  
  Xbar <- mean(X) # Sample mean of the simulated values X
  demeaned_X <- X - Xbar # X - Xbar
  var = mean((X - Xbar)^2) 
  # Berechnet die Differenzen zwischen den beiden Punkten t1 und t2
  lag <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  
  # Berechnet die euklidischen Distanzen zwischen allen Punkten im Raster
  dist_matrix <- as.matrix(dist(grid))
  
  # Anwenden der Kernel-Funktion auf die Distanzen
  K_vals <- kernel_fun(lag - dist_matrix, bandwidth)
  
  X_ij <- outer(demeaned_X, demeaned_X, "*") 
  
  cov_vals = X_ij / var
  
  # Zähler
  numerator <- sum(K_vals * cov_vals)
  
  # Nenner
  denominator <- sum(K_vals)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  return(list(covariance = numerator / denominator, weights = K_vals))
}

#-------------------------------------------------------------------------------
# Kernel-Kovarianzfunktion mit difference

kernel_cov_flo <- function(t1, t2, X, grid, bandwidth,
                           kernel_function="gaussian_kernel") {
  # Ensure the kernel_function is correctly specified
  kernel_fun = switch(
    kernel_function,
    gaussian_kernel = gaussian_kernel,
    epanechnikov_kernel = epanechnikov_kernel,
    rechteck_kernel = rechteck_kernel,
    triangle_kernel = triangle_kernel,
    cubic_b_spline_kernel = cubic_b_spline_kernel,
    stop("Specify either gaussian_kernel or epanechnikov_kernel
         or rechteck_kernel or triangle_kernel kernel 
         or or cubic_b_spline_kernel")
  )
  
  
  # Sample mean of the simulated values X
  Xbar <- mean(X)
  demeaned_X <- X - Xbar # Demeaned X values
  var = mean((X - Xbar)^2) 
  # Number of points in the grid
  n <- length(grid$x)
  
  # Initializing numerator and denominator for covariance calculation
  numerator <- 0
  denominator <- 0
  
  # Looping through all pairs of points (i, j)
  for (i in 1:n) {
    for (j in 1:n) {
      # Difference in coordinates for the grid points
      diff_x <- grid$x[i] - grid$x[j]
      diff_y <- grid$y[i] - grid$y[j]
      
      # Apply the chosen kernel function
      k_x <- kernel_fun(((t1[1] - t2[1]) - diff_x),  bandwidth)
      k_y <- kernel_fun(((t1[2] - t2[2]) - diff_y),  bandwidth)
      K_val <- k_x * k_y
      
      # Covariance components
      numerator <- numerator + K_val * ((demeaned_X[i] * demeaned_X[j]) / var)
      denominator <- denominator + K_val 
    }
  }
  
  # Avoid division by zero
  if (denominator == 0) {
    return(list(covariance = NA, weights = matrix(0, n, n)))
  }
  
  return(list(covariance = numerator / denominator, 
              weights = matrix(K_val, n, n)))
}

#-------------------------------------------------------------------------------

kernel_cov_vec = function(X, grid, bandwidth, 
                          kernel_function = "gaussian_kernel"){
  # Kernel Funktion auswählen
  kernel_fun = switch(
    kernel_function,
    "gaussian_kernel" = gaussian_kernel,
    "epanechnikov_kernel" = epanechnikov_kernel,
    "rechteck_kernel" = rechteck_kernel,
    "triangle_kernel" = triangle_kernel,
    "cubic_b_spline_kernel" = cubic_b_spline_kernel,
    stop("Specify either 'gaussian_kernel' or 'epanechnikov_kernel' or 
         'rechteck_kernel' or 'triangle_kernel' or 'cubic_b_spline_kernel'")
  )
  
  n <- nrow(grid)
  
  Xbar = mean(X)  # Mittelwert von X
  demeaned_X = X - Xbar  # Zentrierte Werte
  var = mean((X - Xbar)^2)  # Varianz der Stichprobe
  
  delta_ij = as.matrix(dist(grid, method = "euclidean", diag = TRUE, upper = TRUE))  # Distanzen
  
  K_vals = matrix(0, n, n)  # Speicher für Kernel-Werte
  cov_vals = matrix(0, n, n)  # Speicher für Kovarianz-Werte
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Kovarianz berechnen
      X_ij <- (X[i] - Xbar) * (X[j] - Xbar)
      cov_vals[i, j] = X_ij / var
      
      # Kernel-Wert berechnen
      K_vals[i, j] = kernel_fun(delta_ij[i, j], bandwidth)
    }
  }
  
  # Zähler und Nenner berechnen
  numerator = sum(K_vals * cov_vals)
  denominator = sum(K_vals)
  
  # Vermeidung von Division durch 0
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  # Ergebnis zurückgeben
  return(list(covariance = numerator / denominator, weights = K_vals))
}

