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
# ANALYTICAL COVARIANCES
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

#-------------------------------------------------------------------------------
# Gneiting covariance function
cov_gneiting <- function(grid, sigma, a, c, beta, method = "euclidean") {
  # COMPUTES THE ANALYTICAL COVARIANCE MATRIX WITH GNEITING
  # Can either be based on eucleadian distance or on the absolute difference
  if(!sum(method == c("difference", "euclidean"))){
    print("Specify either euclidean or difference method")
  }
  
  if(method == "euclidean"){
    # Euclidean distance
    h <- as.matrix(dist(grid, method = "euclidean"))
    
    #covariance = sigma^2 / (((a * abs(h) + 1)^2) * ((b * abs(h))^2 + 1)^1.5)
    
    covariance = (sigma^2 / (a * abs(h) + 1)) * 
      exp(-c * abs(h) / (a * abs(h) + 1)^(beta / 2))
    
  }else if (method == "difference"){
    # Absolute difference
    h1 <- abs(outer(grid$x, grid$x, "-"))
    h2 <- abs(outer(grid$y, grid$y, "-"))
    
    #covariance = sigma^2 / (((a * abs(h2) + 1)^2) * ((b * abs(h1))^2 + 1)^1.5)
    
    covariance = (sigma^2 / (a * h2 + 1)) * 
      exp(-c * h1 / (a * h2 + 1)^(beta / 2))
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
  
  
  # Check that grid has the correct column names "t1" and "t2"
  expected_names <- c("t1", "t2")
  
  if (!all(colnames(grid) == expected_names)) {
    stop("Error: The grid must have column names 't1' and 't2'.")
  }
  
  # Compute pairwise differences for each dimension
  h1_diff = outer(grid$t1, grid$t1, "-")  # Difference in first dimension
  h2_diff = outer(grid$t2, grid$t2, "-")  # Difference in second dimension
  
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
  indicator <- ifelse(abs(x) <= 1, 1, 0)
  
  # Berechne den Kernel-Wert
  kernel_value <- (1 - abs(x)) * indicator / bandwidth
  
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
    stop("Specify either gaussian_kernel or epanechnikov_kernel
         or rechteck_kernel or or triangle_kernel kernel")
  )
  
  Xbar = mean(X) # Sample mean of the simulated values X
  demeaned_X = X - Xbar # X- Xbar
  
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
  n = length(X)
  Xsum =  1/n * sum(X)
  #weight = weight_function(t1,t2)
  
  # Numerator
  numerator = sum(K_vals * (X_ij / Xsum))
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
  kernel_fun <- switch(
    kernel_function,
    gaussian_kernel = gaussian_kernel,
    epanechnikov_kernel = epanechnikov_kernel,
    rechteck_kernel = rechteck_kernel,
    triangle_kernel = triangle_kernel,
    stop("Specify either gaussian_kernel or epanechnikov_kernel
         or rechteck_kernel or or triangle_kernel kernel")
  )
  
  Xbar <- mean(X) # Sample mean of the simulated values X
  demeaned_X <- X - Xbar # X - Xbar
  
  # Berechnet die Differenzen zwischen den beiden Punkten t1 und t2
  lag <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  
  # Berechnet die euklidischen Distanzen zwischen allen Punkten im Raster
  dist_matrix <- as.matrix(dist(grid))
  
  # Anwenden der Kernel-Funktion auf die Distanzen
  K_vals <- kernel_fun(lag - dist_matrix, bandwidth)
  
  X_ij <- outer(demeaned_X, demeaned_X, "*") 
  
  weight = weight_function(t1,t2)
  
  # Zähler
  numerator <- sum(K_vals * X_ij)
  
  # Nenner
  denominator <- sum(K_vals)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  return(list(covariance = numerator * weight / denominator, weights = K_vals))
}

#-------------------------------------------------------------------------------
# Kernel-Kovarianzfunktion mit difference

kernel_cov_flo <- function(t1, t2, X, grid, bandwidth,
                           kernel_function="gaussian_kernel") {
  # Ensure the kernel_function is correctly specified
  kernel_fun <- switch(
    kernel_function,
    gaussian_kernel = gaussian_kernel,
    epanechnikov_kernel = epanechnikov_kernel,
    rechteck_kernel = rechteck_kernel,
    triangle_kernel = triangle_kernel,
    stop("Specify either gaussian_kernel or epanechnikov_kernel
            or rechteck_kernel or triangle_kernel kernel")
  )
  
  # Sample mean of the simulated values X
  Xbar <- mean(X)
  demeaned_X <- X - Xbar # Demeaned X values
  
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
      
      # Apply the weight function
      #weight <- weight_function(t1, t2)
      
      # Covariance components
      numerator <- numerator + K_val * demeaned_X[i] * demeaned_X[j]
      denominator <- denominator + K_val 
    }
  }
  
  # Avoid division by zero
  if (denominator == 0) {
    return(list(covariance = NA, weights = matrix(0, n, n)))
  }
  
  return(list(covariance = numerator / denominator, weights = matrix(K_val, n, n)))
}


#-------------------------------------------------------------------------------
kernel_function <- function(t_diff, h, kernel_type = "triangular") {
  # Kernel function: default is triangular
  if (kernel_type == "triangular") {
    return(pmax(1 - abs(t_diff) / h, 0))
  } else {
    stop(paste("Kernel type", kernel_type, "not recognized"))
  }
}

nonparametric_covariance_estimator <- function(X, t, h, kernel_type = "triangular") {
  # Computes the nonparametric covariance estimator for a stationary random field
  #
  # Args:
  #   X: Matrix of observed data points (n x d) where n is the number of points and d is the dimension
  #   t: Vector of time points where the covariance function is estimated
  #   h: Bandwidth parameter for kernel smoothing
  #   kernel_type: The type of kernel to use
  #
  # Returns:
  #   Vector of estimated covariance values at the time points t
  
  n <- nrow(X)
  d <- ncol(X)
  p_hat <- numeric(length(t))
  
  # Compute the mean of the observed data
  X_mean <- colMeans(X)
  
  # Loop over time points t
  for (k in seq_along(t)) {
    t_k <- t[k]
    numerator_sum <- 0
    denominator_sum <- 0
    
    # Double summation over i and j
    for (i in 1:n) {
      for (j in 1:n) {
        t_diff <- t_k - sqrt(sum((X[i, ] - X[j, ])^2))  # Euclidean distance between points i and j
        kernel_value <- kernel_function(t_diff, h, kernel_type)
        
        # Compute the difference between X[i] and X_mean, and X[j] and X_mean
        X_ij_diff <- (X[i, ] - X_mean) * (X[j, ] - X_mean)
        numerator_sum <- numerator_sum + kernel_value * sum(X_ij_diff)
        denominator_sum <- denominator_sum + kernel_value
      }
    }
    
    # Calculate the estimate p_hat for time t_k
    if (denominator_sum != 0) {
      p_hat[k] <- numerator_sum / denominator_sum
    } else {
      p_hat[k] <- 0  # Handle case where denominator is zero
    }
  }
  
  return(p_hat)
}

# Example usage
set.seed(123)
n <- 100  # Number of observations
d <- 2    # Dimension of the random field
X <- matrix(rnorm(n * d), n, d)  # Simulated data
t <- seq(0, 5, length.out = 100)  # Time points where the covariance is estimated
h <- 0.4  # Bandwidth

p_hat <- nonparametric_covariance_estimator(X, t, h)
print(p_hat)

