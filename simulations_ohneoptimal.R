# MASTERARBEIT DARWIN 
#-------------------------------------------------------------------------------
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)
if (!require(tictoc)) install.packages("tictoc"); library(tictoc)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(pheatmap)) install.packages("pheatmap"); library(pheatmap)
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("utils_functions.R") # Upload functions in separate script
source("plotting.R") # Upload functions in separate script
set.seed(1234)# set seed for replication
#-------------------------------------------------------------------------------
# SIMULATION OF GRID
grid_size = 20
xdim = grid_size # horizontal size
ydim = grid_size #vertical size

# Erzeugt zufällige Koordinaten im Bereich [0, 1] für das Raster
grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))

# ORDERED GRID
#grid = grid %>% arrange(x, y)

# NORM ORDERED GRID
# grid$norm = sqrt(grid$x^2 + grid$y^2)
# grid = grid %>% arrange(norm)
# grid = grid[, 1:2]


# Selecting Points with many neighbours and few neighbours
point_many = select_point_by_neighbour(grid, choice = "many", perc = 0.1)
point_few = select_point_by_neighbour(grid, choice = "few", perc = 0.1)
point_border = select_point_by_neighbour(grid, choice = "border",
                                            perc = 0.2, margin = 0.05)
index_many = point_many$index
index_few = point_few$index
index_border = point_border$index
# #-----------------------------------------------------------------------------
# # GENERATE TRUE ANALYTICAL COVARIANCE
# #-----------------------------------------------------------------------------
# 
sigma = 1
phi = 3
true_covariance = cov_exponential(grid, sigma, phi, method = "difference")
plot_matrix_image(true_covariance, main = paste0("True Covariance Matrix\n",
                                              "phi=", phi,", ",
                                              "sigma=", sigma))
#-------------------------------------------------------------------------------
# # SIMULATE SPATIAL GAUSSIAN RANDOM FIELD
X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_covariance)
grid$sim = X


# Add a 'type' variable to indicate point categories
grid$type = "Normal"
grid$type[index_many] = "Many"
grid$type[index_few] = "Few"
grid$type[index_border] = "Border"

# Define custom shapes for each type
shape_values = c("Normal" = 16, "Many" = 2, "Few" = 3, "Border" = 4)

# Create the plot
ggplot(grid, aes(x = x, y = y, color = sim, shape = type)) +
  geom_point(size = 5) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = shape_values) +
  labs(title = "Simulated Spatial Data",
       x = "X Coordinate",
       y = "Y Coordinate",
       color = "Value",
       shape = "Point Type") +
  theme_minimal()

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

sigma = 1
phi = 2
perc = 0.1
Ns = c(30, 50, 70) # grid sizes
bandwidths = c(0.01, 0.02, 0.03) # bandwidths
kernel = "rechteck_kernel"
truecov_matrices = list()
estcov_matrices = list()

covariances_df = data.frame(
  grid_size = integer(),
  phi = numeric(),
  bandwidth = numeric(),
  spectral_norm = numeric(),
  frob_norm = numeric()
)

point_diff_df = data.frame(
  grid_size = integer(),
  phi = numeric(),
  bandwidth = numeric(),
  point_type = character(),
  distance_type = character(),
  true_cov = numeric(),
  est_cov = numeric(),
  difference = numeric()
)



tic("Loop for different Grid Sizes")
for(N in Ns){
  cat("Grid Size: ", N, "\n")
  estcov_matrices[[as.character(N)]] = list()
  
  # GRID GENERATION #
  xdim = N # horizontal size
  ydim = N #vertical size
  # Erzeugt zufällige Koordinaten im Bereich [0, 1] für das Raster
  grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))
  # ORDERING OF THE GRID (ASCENDING) 
  grid = grid %>% arrange(x, y)
  # Selecting Points with many neighbours and few neighbours
  point_many = select_point_by_neighbour(grid, choice = "many", perc = perc)
  point_few = select_point_by_neighbour(grid, choice = "few", perc = perc)
  point_border = select_point_by_neighbour(grid, choice = "border", perc = 0.1,
                                           margin = 0.05)
  
  points = rbind(point_many$point,
                  point_few$point,
                  point_border$point,
                  make.row.names=FALSE)
  
  # True Covariance
  truecov_matrices[[as.character(N)]] = cov_exponential(grid,
                                                        sigma, phi,
                                                        method = "difference")
  
  # Spatial Simulation
  X = mvrnorm(n = 1, mu = rep(0, nrow(grid)),
              Sigma = truecov_matrices[[as.character(N)]])
  grid$sim = X
  
  tic("Loop for different Bandwidths")
  for(bw in bandwidths){
    cat("Grid Size:", N, "Bandwidths: ", bw, "\n")
    
    est_cov_matrix = matrix(NA, nrow = nrow(grid),
                            ncol = nrow(grid))
    for(i in 1:nrow(grid)) {
      for(j in 1:nrow(grid)) {
        est_cov_matrix[i,j] = kernel_cov(as.numeric(grid[i, 1:2]),
                                         as.numeric(grid[j, 1:2]),
                                         grid$sim, grid[,1:2],
                                         bw, kernel)$covariance
      }
    }
    estcov_matrices[[as.character(N)]][[as.character(bw)]] = est_cov_matrix

  # COVARIANCE
  spectral_norm = norm(truecov_matrices[[as.character(N)]] - 
                         estcov_matrices[[as.character(N)]][[as.character(bw)]],
                       type = "2")
  frob_norm = norm(truecov_matrices[[as.character(N)]] - 
                         estcov_matrices[[as.character(N)]][[as.character(bw)]],
                       type = "F")
  
  covariances_df = rbind(covariances_df, data.frame(
    grid_size = N,
    phi = phi,
    bandwidth = bw,
    spectral_norm = spectral_norm,
    frob_norm = frob_norm
  ))
  
  
  # POINTS
  for(i in 1:nrow(points)){
    distance_points = get_distance_points(grid[, c("x", "y")],
                                         points[i, c("x", "y")],
                                         exclude_self = T)
    x_point = points[i, c("x")]
    y_point = points[i, c("y")]
    type = points[i, c("type")]
    
    for(j in 1:nrow(distance_points)){
      distance_type = distance_points[j, "type"]
      true_point_cov = sigma^2*exp(-(x_point - distance_points[j, "x"] / phi)*
                          exp(-(y_point - distance_points[j, "y"]) / phi))
                  
      est_point_cov = kernel_cov(
                        as.numeric(points[i, c("x", "y")]),
                        as.numeric(distance_points[j, c("x", "y")]),
                        X,
                        grid[,1:2],
                        bw,
                        kernel)
      
      point_diff_df = rbind(point_diff_df, data.frame(
        grid_size = N,
        phi = phi,
        bandwidth = bw,
        point_type = type,
        distance_type = distance_type,
        true_cov = true_point_cov,
        est_cov = est_point_cov$covariance,
        difference = abs(true_point_cov - est_point_cov$covariance)
        )
      )
    }
    }
  }
  toc()
}
toc()
print(covariances_df)
print(point_diff_df)
#-------------------------------------------------------------------------------
# Save the results
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))

timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
file_name_cov = paste0("results_cov_", timestamp, ".csv")
file_name_points = paste0("results_poi_", timestamp, ".csv")

write.csv(covariances_df, file = file.path(wd, file_name_cov),
          row.names = F)
write.csv(covariances_df, file = file.path(wd, file_name_points),
          row.names = F)
