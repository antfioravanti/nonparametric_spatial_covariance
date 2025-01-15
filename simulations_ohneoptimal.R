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

xdim = 50 # horizontal size
ydim = 50 #vertical size

# Erzeugt zufällige Koordinaten im Bereich [0, 1] für das Raster
grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))
grid = grid %>% arrange(x, y)
# Selecting Points with many neighbours and few neighbours
point_many = select_point_by_neighbour(grid, choice = "many", perc = 0.1)
point_few = select_point_by_neighbour(grid, choice = "few", perc = 0.1)
point_border = select_point_by_neighbour(grid, choice = "border",
                                            perc = 0.2, margin = 0.05)
index_many = point_many$index
index_few = point_few$index
index_border = point_border$index
#-------------------------------------------------------------------------------
# GENERATE TRUE ANALYTICAL COVARIANCE
#-------------------------------------------------------------------------------

sigma = 1
phi = 3
true_covariance = cov_exponential(grid, sigma, phi, method = "difference")
plot_matrix_image(true_covariance, main = paste0("True Covariance Matrix\n",
                                              "phi=", phi,", ",
                                              "sigma=", sigma))
#-------------------------------------------------------------------------------
# SIMULATE SPATIAL GAUSSIAN RANDOM FIELD
X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_covariance)
grid$sim = X


# Add a 'type' variable to indicate point categories
grid$type <- "Normal"
grid$type[index_many] <- "Many"
grid$type[index_few] <- "Few"
grid$type[index_border] <- "Border"

# Define custom shapes for each type
shape_values <- c("Normal" = 16, "Many" = 2, "Few" = 3, "Border" = 4)

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

