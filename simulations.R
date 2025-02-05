# MASTERARBEIT DARWIN 
#-------------------------------------------------------------------------------
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)
if (!require(tictoc)) install.packages("tictoc"); library(tictoc)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(pheatmap)) install.packages("pheatmap"); library(pheatmap)
if (!require(ggpubr)) install.packages("ggpubr"); library(ggpubr)
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("utils_functions.R") # Upload functions in separate script
set.seed(1234)# set seed for replication
#-------------------------------------------------------------------------------
# SIMULATION OF GRID

xdim = 100 # horizontal size
ydim = 100 #vertical size

# Erzeugt zufällige Koordinaten im Bereich [0, 1] für das Raster
grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))
#-------------------------------------------------------------------------------
plot1 <- plot(grid$x, grid$y,
     main = "Räumliche Positionen aus einer gleichmäßigen Zufallsverteilung",
     xlab = "x-Koordinate", ylab = "y-Koordinate",
     pch = 16, col = "blue")
text(grid$x, grid$y, labels = 1:nrow(grid), pos = 3, cex = 0.7, col = "red")
#-------------------------------------------------------------------------------
# GENERATE TRUE ANALYTICAL COVARIANCE
#-------------------------------------------------------------------------------



sigma = 1
phi = 1/3
#phi_x = 0.5
#phi_y = 1

#true_covariance = cov_exponential(grid, sigma, phi = NULL, phi_x, phi_y, method = "difference")

true_covariance = cov_exponential(grid, sigma, phi, method = "difference")

titel_cov = paste("True Exponential Covariance;", "sigma =", sigma, "phi =", phi)

# Plotten der Kovarianzmatrix
plot_matrix(true_covariance, main = titel_cov, labels = T)


true_cov_melted <- melt(true_covariance)

# Heatmap der Kovarianzmatrix
ggplot(true_cov_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "blue", mid = "white", midpoint = 0) +
  labs(title = "Heatmap der Kovarianzmatrix", x = "Punktindex", y = "Punktindex") +
  theme_minimal()


# sigma = 1
# alpha1 = 1
# alpha2 = 1
# lambda1 = 5
# lambda2 = 5
# beta = 0
# true_covariance1 = ModifiedExponentialCovariance(grid,sigma,alpha1,alpha2,lambda1
#  ,lambda2, beta, test_sep =T)
# 
# true_cov = true_covariance1$covariance
# 
# # Plotten der Kovarianzmatrix
# plot_matrix(true_cov, main = "Exponential Covariance", labels = F)


#-------------------------------------------------------------------------------

#Check covariance is positive semidefinite
is_positive_semi_definite(true_covariance)

#-------------------------------------------------------------------------------
#SIMULATE SPATIAL GAUSSIAN RANDOM FIELD
X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_covariance)
grid$sim = X

#Plot the process where the colors indicate the value of the spatial process
ggplot(grid, aes(x = x, y = y, color = sim)) +
  geom_point(size = 5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Simuliertes räumliches gaußsches Zufallsfeld", x = "x-Koordinate", y = "y-Koordinate",
       color = "Wert") +
  theme_minimal()



#X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_covariance)

#-------------------------------------------------------------------------------
# ESTIMATE SPATIAL COVARIANCE

bandwidth = 0.01 # Beispiel-Bandbreite

#kernel = "gaussian_kernel"

#kernel = "epanechnikov_kernel"

kernel = "rechteck_kernel"

#kernel = "triangle_kernel"

#kernel = "cubic_b_spline_kernel"

#-------------------------------------------------------------------------------

tic("matrix cov estimation")
est_cov_matrix = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
for(i in 1:nrow(grid)) {
  for(j in 1:nrow(grid)) {
    est_cov_matrix[i,j] = kernel_cov(as.numeric(grid[i, 1:2]),
                                     as.numeric(grid[j, 1:2]),
                                     X, grid[,1:2],
                                     bandwidth, kernel)$covariance
  }
}
toc()
 
is_positive_semi_definite(est_cov_matrix)

norm_value <- norm(est_cov_matrix - true_covariance, type = "2")
print(paste("NORM: ", norm_value))


mse_error <- mean((est_cov_matrix - true_covariance)^2)
print(paste("MSE Error: ", mse_error))


# tic("vector cov estimation")
# est_cov_vector = numeric(nrow(grid))
# 
# # Schleife über jeden Punkt im Gitter
# for (i in 1:nrow(grid)) {
#   # Schätze die Kovarianz für den Punkt grid[i, 1:2]
#   est_cov = kernel_cov_new(as.numeric(grid[i, 1:2]), X, grid[, 1:2], bandwidth,
#                            kernel)
#   est_cov_vector[i] = est_cov$covariance
# }
# toc()
# View(est_cov_vector)


#-------------------------------------------------------------------------------

# Prüfen für nur ein Punktpaar auf dem Gitter

result <- kernel_cov(as.numeric(grid[22, 1:2]),
           as.numeric(grid[33, 1:2]),
           X, grid[,1:2],
           bandwidth, kernel)

result1 <- kernel_cov(as.numeric(grid[101, 1:2]),
                     as.numeric(grid[110, 1:2]),
                     X, grid[,1:2],
                     bandwidth, kernel)

result2 <- kernel_cov(as.numeric(grid[72, 1:2]),
                     as.numeric(grid[64, 1:2]),
                     X, grid[,1:2],
                     bandwidth, kernel)

result3 <- kernel_cov(as.numeric(grid[1, 1:2]),
                      as.numeric(grid[1, 1:2]),
                      X, grid[,1:2],
                      bandwidth, kernel)


View(result$weights)

#-------------------------------------------------------------------------------

threshold <- 0.1

# Créer une copie de la matrice originale pour la nouvelle version
est_cov_matrix_new <- est_cov_matrix
true_covariance_new <- true_covariance

# Appliquer le seuil uniquement hors diagonale
est_cov_matrix_new[abs(est_cov_matrix_new) > threshold & row(est_cov_matrix_new)
                   != col(est_cov_matrix_new)] <- 0

#true_covariance_new[abs(true_covariance_new) > threshold & row(true_covariance_new)
                #!= col(true_covariance_new)] <- 0

# Calculer la norme pour vérifier la convergence
norm_value <- norm(est_cov_matrix_new - true_covariance_new, type = "2")
print(norm_value)

#-------------------------------------------------------------------------------

threshold <- 0.1

# Appliquer le seuil sur les valeurs absolues
est_cov_matrix_new[abs(est_cov_matrix_new) > threshold] <- 0
true_covariance_new[abs(true_covariance_new) > threshold] <- 0

# Calculer la norme
norm_value <- norm(est_cov_matrix_new - true_covariance_new, type = "2")
print(norm_value)

#-------------------------------------------------------------------------------

threshold <- 0.03

# Set large values in both matrices to zero
est_cov_matrix_new[est_cov_matrix_new > threshold] <- 0
#true_covariance_new[true_covariance_new > threshold] <- 0

# Compute the norm to check convergence
norm_value <- norm(est_cov_matrix_new - true_covariance_new, type = "2")
print(norm_value)

#-------------------------------------------------------------------------------

# tic(" auto matrix cov estimation") # Berechne die Autokorrelationsmatrix
# est_corr_matrix = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
# for(i in 1:nrow(grid)) {
#   for(j in 1:nrow(grid)) {
#     est_corr_matrix[i,j] = est_cov_matrix[i,j] /
#       sqrt(est_cov_matrix[i,i] * est_cov_matrix[j,j])
#   }
# }
# toc()
# 
# 
# # Prüfen, ob die Diagonalelemente 1 sind
# diag(est_corr_matrix)
# is_positive_semi_definite(est_corr_matrix)


#-------------------------------------------------------------------------------
# Norm 12.17411
norm(est_cov_matrix - true_covariance, type = "2")

#-------------------------------------------------------------------------------
# Norm with autokorellation
norm(est_corr_matrix - true_covariance, type = "2")


tic("matrix cov estimation")
est_cov_matrix_flo = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
for(i in 1:nrow(grid)) {
  for(j in 1:nrow(grid)) {
    est_cov_matrix_flo[i,j] = kernel_cov_flo(as.numeric(grid[i, 1:2]),
                                     as.numeric(grid[j, 1:2]),
                                     X, grid[,1:2],
                                     bandwidth, kernel)$covariance
  }
}
toc()


is_positive_semi_definite(est_cov_matrix_flo)

norm(est_cov_matrix_flo - true_covariance, type = "2")

#-------------------------------------------------------------------------------

# Differenzmatrix berechnen
diff_matrix <- est_cov_matrix - true_covariance

# Heatmap der Differenzmatrix erstellen
pheatmap(diff_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap der Differenz")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


grid_sizes <- seq(10, 100, by = 10) # Verschiedene Gittergrößen
mse_values <- numeric(length(grid_sizes)) # MSE-Werte speichern

for (i in seq_along(grid_sizes)) {
  xdim <- grid_sizes[i]
  ydim <- grid_sizes[i]
  grid <- data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))
  
  # Erzeuge die wahre Kovarianzmatrix (oder eine Schätzungsfunktion)
  true_covariance <- cov_exponential(grid, sigma, phi, method = "difference")
  X = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_covariance)
  grid$sim = X
  
  # Schätze die Kovarianz mit deinem Kernel
  est_cov_matrix <- matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
  for(i in 1:nrow(grid)) {
    for(j in 1:nrow(grid)) {
      est_cov_matrix[i,j] = kernel_cov(as.numeric(grid[i, 1:2]),
                                       as.numeric(grid[j, 1:2]),
                                       X, grid[,1:2],
                                       bandwidth, kernel)$covariance
    }
  }
  
  # Berechne den MSE
  mse_values[i] <- mean((est_cov_matrix - true_covariance)^2)
}

# MSE visualisieren
plot(grid_sizes, mse_values, type = "b", main = "MSE vs Grid Size",
     xlab = "Grid Size", ylab = "MSE", pch = 19, col = "blue")

