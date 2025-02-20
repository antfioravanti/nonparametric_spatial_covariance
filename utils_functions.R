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
  #   alpha1 and alpha2 control the smoothness in each dimension
  #   lambda1 and lamda2 control the range in each dimension
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
    cov_sep = sigma^2 * exp(-(abs(h1_diff)^alpha1 / lambda1)-beta*abs(h1_diff))* 
      exp(-(abs(h2_diff)^alpha2 / lambda2) -beta*abs(h2_diff))
    
    norm_diff = norm(cov_sep - cov, type = "2")
    return(list(covariance = cov, norm_diff = norm_diff))
  }else{
    return(cov)
  }
}


#-------------------------------------------------------------------------------
# GRID
#-------------------------------------------------------------------------------

select_point_by_neighbour = function(grid,
                                     choice = c("many", "few", "border"),
                                     perc = 0.2,
                                     margin = 0.05,
                                     randomize = TRUE) {
  # grid: A data.frame with 'x' and 'y' columns representing points
  # choice: 
  #   "many" to select a point with many neighbors,
  #   "few" for few neighbors,
  #   "border" for an isolated point near the border
  # perc: Percentage of Maximum Distance to set threshold
  # margin: Proximity to border (0 < margin < 0.5) for "border" choice
  # randomize: If TRUE, randomly select among multiple candidates
  
  # Validate input - MAKE SURE YOU HAVE x and y IN THE COLUMNS OF THE GRID
  if (!all(c("x", "y") %in% names(grid))) {
    stop("Grid must have 'x' and 'y' columns.")
  }
  
  choice = match.arg(choice)
  
  # Pairwise distances - euclidean distance is fine here
  distance_matrix = as.matrix(dist(grid[, c("x", "y")]))
  
  # Threshold distance
  threshold = perc * max(distance_matrix)
  
  if (choice == "border") {
    # Identify border points based on margin - we allow for a margin of 0.05 
    is_border = with(grid, x <= margin | x >= (1 - margin) | 
                       y <= margin | y >= (1 - margin))
    border_points = grid[is_border, ]
    
    if (nrow(border_points) == 0) {
      stop("No points found near the border with the specified margin.")
    }
    
    # Pairwise distances among border points
    distance_matrix_border = as.matrix(dist(border_points[, c("x", "y")]))
    
    # Threshold distance for border points
    threshold_border = perc * max(distance_matrix_border)
    
    # Count neighbors within threshold for each border point (excluding itself)
    neighbor_counts_border = apply(distance_matrix_border, 1, 

    function(row) sum(row < threshold_border) - 1)

    
    # Counting border_points
    border_points$neighbor_count = neighbor_counts_border
    
    # Find minimum neighbor count among border points
    min_neighbors = min(border_points$neighbor_count)
    
    # Identify border points with minimum neighbor count
    isolated_border_points = border_points[border_points$neighbor_count == min_neighbors, ]
    
    if (nrow(isolated_border_points) == 0) {
      stop("No isolated border points found with the specified threshold.")
    }
    
    # Select one point: 
    if (nrow(isolated_border_points) > 1 && randomize) {
      # randomly
      selected_row = isolated_border_points[sample(nrow(isolated_border_points), 1), ]
    } else {
      # first point in the list
      selected_row = isolated_border_points[1, ]
    }
    
    # Find the index of the selected point in the original grid
    selected_index = which(grid$x == selected_row$x & grid$y == selected_row$y)
    selected_row$type = "Border"
    # Return the selected point and its index
    return(list(point = selected_row[,c("x", "y", "type")], 
                index = selected_index))
    
  } else {
    # For "many" and "few" choices
    
    # Count neighbors within threshold for each point (excluding itself)
    neighbor_counts = apply(distance_matrix, 1, 
                            function(row) sum(row < threshold) - 1)
    
    if (choice == "many") {
      selected_index = which.max(neighbor_counts)
      grid$type = "Many"
    } else { # choice == "few"
      selected_index = which.min(neighbor_counts)
      grid$type = "Few"
    }
    
    # Return the selected point and its index
    return(list(point = grid[selected_index, ], index = selected_index))
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
  return(ifelse(abs(x / bandwidth) <= 1, (1 - abs(x / bandwidth)), 0))
}

cubic_b_spline_kernel <- function(x, bandwidth) {
  # CUBIC B-SPLINE KERNEL FUNKTION
  u <- x / bandwidth
  kernel_value <- 0.5 * exp(-abs(u) / sqrt(2)) * (sin(abs(u) / sqrt(2)) + pi /4)
  return(kernel_value)
}

#-------------------------------------------------------------------------------
# Estimators 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
weight_function <- function(t1,t2){
  log <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  ifelse(abs(log) <= 10, 1, 0)
}
#-------------------------------------------------------------------------------

# Vektorisierte Kernel-Kovarianzfunktion mit Differenzen

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
  #sample_var = mean((X - Xbar)^2)  # Varianz der Stichprobe
  
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
  
  #cov_vals = X_ij / sample_var
  
  
  # Numerator
  numerator = sum(K_vals * X_ij)
  # Denominator
  denominator = sum(K_vals)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  # Calcul du nombre de 0 et 1 dans les poids (K_vals)
  count_zeros = sum(K_vals == 0) # Compte les 0 dans la matrice
  count_ones = sum(K_vals == (1 / bandwidth^2)) # Compte les 1 / bandwidth^2 dans la matrice
  
  return(list(covariance = numerator / denominator, 
              weights = K_vals, 
              count_zeros = count_zeros, 
              count_ones = count_ones))
}

get_distance_points <- function(grid, point_df, exclude_self = TRUE) {
  # grid: A data.frame with 'x' and 'y' columns representing points
  # point_df: A data.frame with exactly one row containing 'x' and 'y' columns
  # exclude_self: If TRUE, excludes the reference point from the grid if present
  
  # Input Validation
  if (!all(c("x", "y") %in% names(grid))) {
    stop("The grid must contain 'x' and 'y' columns.")
  }
  
  if (!is.data.frame(point_df) || nrow(point_df) != 1 || 
      !all(c("x", "y") %in% names(point_df))) {
    stop("The reference point must be a data.frame with exactly one row and 'x' and 'y' columns.")
  }
  
  # Extract reference coordinates
  ref_x <- point_df$x[1]
  ref_y <- point_df$y[1]
  
  # Calculate Euclidean distances from the reference point to all grid points
  distances <- sqrt((grid$x - ref_x)^2 + (grid$y - ref_y)^2)
  
  # Optionally exclude the reference point itself from the grid
  if (exclude_self) {
    same_point_indices <- which(grid$x == ref_x & grid$y == ref_y)
    if (length(same_point_indices) > 0) {
      distances[same_point_indices] <- Inf  # Assign infinite distance to exclude
    }
  }
  
  # Identify the Closest Point
  closest_index <- which.min(distances)
  closest_point <- grid[closest_index, ]
  
  # Identify the Medium-Distance Point
  median_distance <- median(distances)
  median_diff <- abs(distances - median_distance)
  medium_index <- which.min(median_diff)
  medium_point <- grid[medium_index, ]
  
  # Identify the Farthest Point
  farthest_index <- which.max(distances)
  farthest_point <- grid[farthest_index, ]
  
  # Compile the results into a single data frame
  result <- data.frame(
    type = c("Closest", "Medium", "Farthest"),
    x = c(closest_point$x, medium_point$x, farthest_point$x),
    y = c(closest_point$y, medium_point$y, farthest_point$y),
    index = c(closest_index, medium_index, farthest_index)
  )
  
  return(result)
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
  #sample_var = mean((X - Xbar)^2) 
  
  # Berechnet die Differenzen zwischen den beiden Punkten t1 und t2
  lag <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  
  # Berechnet die euklidischen Distanzen zwischen allen Punkten im Raster
  #dist_matrix <- as.matrix(dist(grid))
  
  # Anwenden der Kernel-Funktion auf die Distanzen
  #K_vals <- kernel_fun(lag - dist_matrix, bandwidth)
  
  K_vals <- kernel_fun(lag, bandwidth)
  
  X_ij <- outer(demeaned_X, demeaned_X, "*") 
  
  #cov_vals = X_ij / sample_var
  
  # Zähler
  numerator <- sum(K_vals * X_ij)
  
  # Nenner
  denominator <- sum(K_vals)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  return(list(covariance = numerator / denominator, weights = K_vals))
}

#-------------------------------------------------------------------------------

# Nicht Vektorisierte Kernel-Kovarianzfunktion mit Differenzen

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
  #sample_var = mean((X - Xbar)^2) 
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
      
      #diff_xy <- sqrt(diff_x^2 + diff_y^2)
      
      # Apply the chosen kernel function
      k_x <- kernel_fun(((t1[1] - t2[1]) - diff_x),  bandwidth)
      k_y <- kernel_fun(((t1[2] - t2[2]) - diff_y),  bandwidth)
      K_vals <- k_x * k_y
      
      # Covariance components
      numerator <- numerator + K_vals * 
        ((demeaned_X[i] * demeaned_X[j]))
      denominator <- denominator + K_vals 
    }
  }
  
  # Avoid division by zero
  if (denominator == 0) {
    return(list(covariance = NA, weights = matrix(0, n, n)))
  }
  
  # Calcul du nombre de 0 et 1 dans les poids (K_vals)
  count_zeros = sum(K_vals == 0) # Compte les 0 dans la matrice
  count_ones = sum(K_vals == (1 / bandwidth^2)) 
  
  return(list(covariance = numerator / denominator, 
              weights = matrix(K_vals, n, n), 
              count_zeros = count_zeros, 
              count_ones = count_ones))
}

#-------------------------------------------------------------------------------

# Nicht Vektorisierte Kernel-Kovarianzfunktion mit euklidischen Distanz

kernel_cov_flo_eucli <- function(t1, t2, X, grid, bandwidth,
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
  sample_var = mean((X - Xbar)^2) 
  # Number of points in the grid
  n <- length(grid$x)
  
  # Initializing numerator and denominator for covariance calculation
  numerator <- 0
  denominator <- 0
  
  lag <- sqrt((t1[1] - t2[1])^2 + (t1[2] - t2[2])^2)
  
  # Looping through all pairs of points (i, j)
  for (i in 1:n) {
    for (j in 1:n) {
      
      K_val <- kernel_fun(lag, bandwidth)
      
      # Covariance components
      numerator <- numerator + K_val * 
        ((demeaned_X[i] * demeaned_X[j]) / sample_var)
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

# Nicht Vektorisierte Kernel-Kovarianzfunktion (Andere Schreibweise ??????)

kernel_cov_vec = function(X, grid, bandwidth, 
                          kernel_function = "gaussian_kernel") {
  
  # Kernel function selection
  kernel_fun = switch(
    kernel_function,
    "gaussian_kernel" = gaussian_kernel,
    "epanechnikov_kernel" = epanechnikov_kernel,
    "rechteck_kernel" = rechteck_kernel,
    "triangle_kernel" = triangle_kernel,
    "cubic_b_spline_kernel" = cubic_b_spline_kernel,
    stop("Specify either 'gaussian_kernel' or 'epanechnikov_kernel' 
          or 'rechteck_kernel' or 'triangle_kernel' or 'cubic_b_spline_kernel'")
  )
  
  n <- nrow(grid)
  Xbar = mean(X)  # Mean of X
  demeaned_X = X - Xbar  # Centered values
  sample_var = mean((demeaned_X)^2)  # Sample variance
  
  K_vals = matrix(0, n, n)  # To store kernel values
  cov_vals = matrix(0, n, n)  # To store covariance values
  
  # Loop through the points in the grid
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate distance between the current pair of grid points
      delta_ij = sqrt(sum((grid[i, 1:2] - grid[j, 1:2])^2))
      
      # Covariance between centered values of X
      X_ij <- (demeaned_X[i]) * (demeaned_X[j])
      cov_vals[i, j] = X_ij / sample_var
      
      # Kernel function value based on the distance
      K_vals[i, j] = kernel_fun(delta_ij, bandwidth)
    }
  }
  
  # Calculate numerator and denominator
  numerator = sum(K_vals * cov_vals)
  denominator = sum(K_vals)
  
  # Avoid division by 0
  if (denominator == 0) return(list(covariance = NA, weights = K_vals)) 
  
  # Return covariance estimate and weights
  return(list(covariance = numerator / denominator, weights = K_vals))
}

#-------------------------------------------------------------------------------

# Vektorisierte Kernel-Kovarianzfunktion (Andere Schreibweise ??????) 

kernel_cov_vec_vectorized <- function(grid, X, bandwidth, kernel_function = "gaussian_kernel") {
  
  # Kernel function selection
  kernel_fun <- switch(
    kernel_function,
    "gaussian_kernel" = gaussian_kernel,
    "epanechnikov_kernel" = epanechnikov_kernel,
    "rechteck_kernel" = rechteck_kernel,
    "triangle_kernel" = triangle_kernel,
    "cubic_b_spline_kernel" = cubic_b_spline_kernel,
    stop("Specify either 'gaussian_kernel', 'epanechnikov_kernel',
         'rechteck_kernel', 'triangle_kernel', or 'cubic_b_spline_kernel'")
  )
  
  # Mean and variance of X
  Xbar <- mean(X)  
  demeaned_X <- X - Xbar  # Centered values
  sample_var <- mean((demeaned_X)^2)  # Sample variance
  
  n <- nrow(grid)
  
  # Compute pairwise distances between all grid points (vectorized)
  distances <- as.matrix(dist(grid[, 1:2]))  # Euclidean distances between all pairs
  
  # Compute pairwise covariances for all pairs (vectorized)
  cov_vals <- outer(demeaned_X, demeaned_X, "*") / sample_var
  
  # Apply the kernel function to the distances (vectorized)
  K_vals <- kernel_fun(distances, bandwidth)
  
  # Compute numerator and denominator for covariance estimation
  numerator <- sum(K_vals * cov_vals)
  denominator <- sum(K_vals)
  
  # Avoid division by 0
  if (denominator == 0) return(list(covariance = NA, weights = K_vals))
  
  # Return the covariance estimate and kernel weights
  return(list(covariance = numerator / denominator, weights = K_vals))
}

#-------------------------------------------------------------------------------

# Fonction pour afficher la covariance pour différentes valeurs de phi
visualize_covariance = function(grid, sigma, phi_values, method = "euclidean") {
  par(mfrow = c(1, length(phi_values)))  # Pour afficher plusieurs graphes côte à côte
  
  for (phi in phi_values) {
    # Calcul de la covariance pour un phi donné
    covariance_matrix = cov_exponential(grid, sigma, phi, method = method)
    
    # Affichage de la carte de chaleur
    image(1:nrow(covariance_matrix), 1:ncol(covariance_matrix), covariance_matrix, 
          main = paste("Covariance (phi =", round(phi, 2), ")"), 
          xlab = "Index x", ylab = "Index y", 
          col = heat.colors(100))
  }
}
