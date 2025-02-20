# MASTERARBEIT DARWIN 
#-------------------------------------------------------------------------------
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)
if (!require(doSNOW)) install.packages("doSNOW"); library(doSNOW)
if (!require(tictoc)) install.packages("tictoc"); library(tictoc)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(pheatmap)) install.packages("pheatmap"); library(pheatmap)
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("utils_functions.R") # Upload functions in separate script
source("plotting.R") # Upload functions in separate script
set.seed(1234)# set seed for replication
#-------------------------------------------------------------------------------
# SETUP OF PARALLEL COMPUTATION

# Create a cluster
num_cores = detectCores() - 5
cl = makeCluster(num_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  source("utils_functions.R") 
  source("plotting.R")
})

#-------------------------------------------------------------------------------
# PARAMETERS
Ns = c(30, 50, 70)               # Grid sizes

#True Covariance
sigma = 1
phi = 2

# Things for point selections
perc = 0.1             # Percentage for selecting many and few neighbors
margin_border = 0.05             # Margin for border selection

# Kernels
bandwidths = c(0.1, 0.2, 0.3)    # Bandwidths
kernel_type = "rechteck_kernel"  # Kernel type

# Initialize a progress bar
pb = txtProgressBar(max = length(bandwidths) * length(Ns), style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

# Initialize data frames to store results
spectral_norm_df = data.frame(
  grid_size = integer(),
  phi = numeric(),
  bandwidth = numeric(),
  spectral_norm = numeric(),
  stringsAsFactors = FALSE
)

cov_diff_df = data.frame(
  grid_size = integer(),
  phi = numeric(),
  bandwidth = numeric(),
  point_type = character(),
  covariance_difference = numeric(),
  stringsAsFactors = FALSE
)

#-------------------------------------------------------------------------------
# Start the parallel simulation
simulation_results = foreach(bw = bandwidths,
                             .combine = rbind,
                             .packages = c("MASS", "dplyr"),
                             .options.snow = opts) %dopar% {
   
   # Initialize temporary data frames for this bandwidth
   temp_spectral_norm = data.frame(
     grid_size = integer(),
     bandwidth = numeric(),
     spectral_norm = numeric(),
     stringsAsFactors = FALSE
   )
   
   temp_cov_diff = data.frame(
     grid_size = integer(),
     bandwidth = numeric(),
     point_type = character(),
     covariance_difference = numeric(),
     stringsAsFactors = FALSE
   )
   
   # Iterate over different grid sizes
   for(N in Ns){
     
     
     # GRID GENERATION #
     xdim = N  # horizontal size
     ydim = N  # vertical size
     
     # Generate random uniformly distributed points in [0, 1] for both x and y
     grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))
     
     # ORDERING OF THE GRID (ASCENDING) 
     grid = grid %>% arrange(x, y)
     
     # Selecting Points with many neighbours, few neighbours, and border
     point_many = select_point_by_neighbour(grid, choice = "many", perc = perc)
     point_few = select_point_by_neighbour(grid, choice = "few", perc = perc)
     point_border = select_point_by_neighbour(grid, choice = "border", 
                                              perc = perc,
                                              margin = margin_border)
     
     # Extract indices
     index_many = point_many$index
     index_few = point_few$index
     index_border = point_border$index
     
     # True Covariance
     true_cov = cov_exponential(grid, sigma, phi, method = "difference")
     
     # Spatial Simulation
     X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_cov)
     grid$sim = X
     
     # Estimate Covariance Matrix using kernel_cov
     # Assuming kernel_cov is defined in utils_functions.R and returns a list with 'covariance'
     est_cov_matrix = matrix(NA, nrow = nrow(grid),
                             ncol = nrow(grid))
     for(i in 1:nrow(grid)) {
       for(j in 1:nrow(grid)) {
         est_cov_matrix[i,j] = kernel_cov(as.numeric(grid[i, 1:2]),
                                          as.numeric(grid[j, 1:2]),
                                          grid$sim, grid[,1:2],
                                          bw, kernel_type)$covariance
       }
     }
     
     # Compute spectral norm difference
     spectral_norm_val = norm(true_cov - est_cov_matrix, type = "2")
     
     # Append to spectral norm data frame
     temp_spectral_norm = rbind(temp_spectral_norm, data.frame(
       grid_size = N,
       bandwidth = bw,
       spectral_norm = spectral_norm_val,
       stringsAsFactors = FALSE
     ))
     
     # IDENTIFY CLOSEST, MEDIUM, AND FAR POINTS FOR SELECTED POINTS
     
     # Combine selected points into a data frame
     selected_points_df = grid[c(index_many, index_few, index_border), 
                               c("x", "y")]
     selected_types = c("Many", "Few", "Border")
     
     # Iterate over each selected point
     for(k in 1:nrow(selected_points_df)) {
       selected_point = selected_points_df[k, ]
       selected_type = selected_types[k]
       
       # Get Closest, Medium, Farthest points relative to the selected point
       reference_point_df = selected_point
       distance_points = get_distance_points(grid, reference_point_df,
                                             exclude_self = TRUE)
       
       # Iterate over each distance-based point
       for(l in 1:nrow(distance_points)) {
         distance_type = distance_points$type[l]
         distance_point = distance_points[l, c("x", "y")]
         
         # Find index_distance by matching coordinates
         index_distance = which(grid$x == distance_point$x & grid$y == distance_point$y)
         
         # Handle cases with multiple matches or no matches
         if(length(index_distance) == 0){
           warning(paste("No matching point found for coordinates:", distance_point$x, distance_point$y))
           next
         }
         if(length(index_distance) > 1){
           index_distance = index_distance[1]  # Take the first match
         }
         
         # Find index_selected by matching coordinates
         index_selected = which(grid$x == selected_point$x & grid$y == selected_point$y)
         
         # Compute covariance difference
         cov_true = true_cov[index_distance, index_selected]
         cov_est = est_cov_matrix[index_distance, index_selected]
         cov_diff = abs(cov_true - cov_est)
         
         # Append to covariance difference data frame
         temp_cov_diff = rbind(temp_cov_diff, data.frame(
           grid_size = N,
           bandwidth = bw,
           point_type = paste(selected_type, distance_type, sep = "_"),
           covariance_difference = cov_diff,
           stringsAsFactors = FALSE
         ))
       }
     }
     
   }
   
   # Combine temporary data frames into a list
   return(list(
     spectral_norm = temp_spectral_norm,
     cov_diff = temp_cov_diff
   ))
}

# Close the progress bar
close(pb)

# Stop the cluster
stopCluster(cl)

#-------------------------------------------------------------------------------
# EXTRACTING RESULTS FROM simulation_results

# Iterate through simulation_results to populate the data frames
for(res in simulation_results){
  spectral_norm_df = rbind(spectral_norm_df, res$spectral_norm)
  cov_diff_df = rbind(cov_diff_df, res$cov_diff)
}

#-------------------------------------------------------------------------------
# VISUALIZATION

# Spectral Norm Differences Plot
ggplot(spectral_norm_df, aes(x = factor(grid_size), y = spectral_norm, fill = factor(bandwidth))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Spectral Norm Differences",
       x = "Grid Size",
       y = "Spectral Norm Difference",
       fill = "Bandwidth") +
  theme_minimal()

# Covariance Differences Plot
ggplot(cov_diff_df, aes(x = point_type, y = covariance_difference, fill = point_type)) +
  geom_boxplot() +
  labs(title = "Covariance Differences for Selected Points",
       x = "Point Type",
       y = "Covariance Difference",
       fill = "Point Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
