### 1 DIMENSIONAL SIMULATION
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("utils_functions.R") # Upload functions in separate script
set.seed(42)
# 1-dimensional simulation as in the paper from Peter Hall and Patil
#-------------------------------------------------------------------------------
kernel_cov_1d = function(t1, t2, X, grid, bandwidth,
                         kernel_function="gaussian_kernel"){
  # COMPUTES THE COVARIANCE ESTIMATOR FOR 1 DIMENSIONAL PROCESS
  
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
  lag = t1 - t2
  # Berechnet die Differenzen zwischen allen x-Koordinaten der Punkte im Raster
  diff_t = outer(grid$t, grid$t, "-")
  
  # Apply the chosen kernel
  k_t = kernel_fun(lag - diff_t, bandwidth)
  
  X_ij = outer(demeaned_X, demeaned_X, "*") 
  
  #weight = weight_function(t1,t2)
  
  # Numerator
  numerator = sum(k_t * X_ij)
  # Denominator
  denominator = sum(k_t)
  
  # Verhindert eine Division durch Null, falls der Nenner Null ist
  if (denominator == 0) return(list(covariance = NA, weights = k_t)) 
  
  return(list(covariance = numerator / denominator, weights = k_t))
  
}

exp_cov = function(grid){
  # EXPONENTIAL COVARIANCE AS SPECIFIED IN THE PAPER
  return(exp(-as.matrix(dist(grid, method = "euclidean")) * 1.98))
}
#-------------------------------------------------------------------------------
tdim = 200
grid = data.frame(t = runif(tdim, 0, 40))
bandwidth = 0.4 # as specified in the paper

# True covariance
true_cov_paper = exp_cov(grid)
is_positive_semi_definite(true_cov_paper)

# Calculate the pairwise distances
pairwise_distances = as.matrix(dist(grid, method = "euclidean"))

# Extract unique distances and corresponding covariances
unique_distances = pairwise_distances[lower.tri(pairwise_distances)]
unique_covariances = true_cov_paper[lower.tri(true_cov_paper)]

X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_cov_paper)
rho_t = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))

for(i in 1:nrow(grid)){
  for(j in 1:nrow(grid)){
    rho_t[i,j] = kernel_cov_1d(as.numeric(grid[i,]), as.numeric(grid[j, ]),
                               X, grid, bandwidth,
                               kernel_function = "triangle_kernel")$covariance
  }}
is_positive_semi_definite(rho_t)


unique_covariances_est = rho_t[lower.tri(rho_t)]

# Create a data frame 
data = data.frame(
  Distance = rep(unique_distances, 2),
  Covariance = c(unique_covariances, unique_covariances_est),
  Type = rep(c("True", "Estimated"), each = length(unique_distances))
)

# Plot true covariance vs estimated covariance as a function of the distance
ggplot(data, aes(x = Distance, y = Covariance, color = Type)) +
  geom_line(size = 0.8) +  # Set the line size to make it thinner
  labs(x = "Distance", y = "Covariance", title = "Covariance vs Distance") +
  theme_minimal() +  # Use a minimal theme for cleaner look
  scale_color_manual(values = c("blue", "red")) +  # Set colors for the lines
  theme(legend.title = element_blank())  # Remove legend title
