################################################################################
################# Erwartungswert Nenner ########################################
################################################################################
library(ggplot2)

# Definieren der zu zeichnenden Funktion
calculate_expression <- function(h1, h2, b, n) {
  
  # Funktion zur Behandlung von Indikatoren mit Bedingung
  indicator <- function(condition) {
    if (condition) return(1) else return(0)
  }
  
  # f1-Funktion mit Indikatorfunktion
  f12 <- ((-0.5 * (h2 + b / 2)^2 + (h2 + b / 2)) * 
            indicator(h2 - b / 2 < 0 & 0 <= h2 + b / 2 & h2 + b / 2 <= 1)) + 
    ((-0.5 * (h2 + b / 2)^2 + (h2 + b / 2) + 
        0.5 * (h2 - b / 2)^2 - (h2 - b / 2)) * 
       indicator(h2 - b / 2 >= 0 & h2 + b / 2 <= 1)) + 
    ((0.5 + 0.5 * (h2 - b / 2)^2 - (h2 - b / 2)) * 
       indicator(0 <= h2 - b / 2 & h2 - b / 2 <= 1 & h2 + b / 2 > 1))
  
  # f2-Funktion mit Indikatorfunktion
  f22 <- ((-0.5 * (h2 - b / 2)^2 - h2 + b / 2) * 
            indicator(-1 <= h2 - b / 2 & h2 - b / 2 < 0 & h2 + b / 2 > 0)) + 
    ((0.5 * (h2 + b / 2)^2 + (h2 + b / 2) - 
        0.5 * (h2 - b / 2)^2 - (h2 - b / 2)) * 
       indicator(h2 - b / 2 >= -1 & h2 + b / 2 <= 0)) + 
    ((0.5 * (h2 + b / 2)^2 + (h2 + b / 2) + 0.5) * 
       indicator(h2 - b / 2 < -1 & -1 <= h2 + b / 2 & h2 + b / 2 <= 0))
  
  # f1-Funktion mit Indikatorfunktion
  f11 <- ((-0.5 * (h1 + b / 2)^2 + (h1 + b / 2)) * 
            indicator(h1 - b / 2 < 0 & 0 <= h1 + b / 2 & h1 + b / 2 <= 1)) + 
    ((-0.5 * (h1 + b / 2)^2 + (h1 + b / 2) + 
        0.5 * (h1 - b / 2)^2 - (h1 - b / 2)) * 
       indicator(h1 - b / 2 >= 0 & h1 + b / 2 <= 1)) + 
    ((0.5 + 0.5 * (h1 - b / 2)^2 - (h1 - b / 2)) * 
       indicator(0 <= h1 - b / 2 & h1 - b / 2 <= 1 & h1 + b / 2 > 1))
  
  # f2-Funktion mit Indikatorfunktion
  f21 <- ((-0.5 * (h1 - b / 2)^2 - h1 + b / 2) * 
            indicator(-1 <= h1 - b / 2 & h1 - b / 2 < 0 & h1 + b / 2 > 0)) + 
    ((0.5 * (h1 + b / 2)^2 + (h1 + b / 2) - 
        0.5 * (h1 - b / 2)^2 - (h1 - b / 2)) * 
       indicator(h1 - b / 2 >= -1 & h1 + b / 2 <= 0)) + 
    ((0.5 * (h1 + b / 2)^2 + (h1 + b / 2) + 0.5) * 
       indicator(h1 - b / 2 < -1 & -1 <= h1 + b / 2 & h1 + b / 2 <= 0))
  
  # Berechnung des Terms
  term0 <- n * indicator((h1 >= -b/2 & h1 <= b/2) | (h2 >= -b/2 & h2 <= b/2))
  #term0 <- n * indicator((h1 >= -b/2 & h1 <= b/2) & (h2 >= -b/2 & h2 <= b/2))
  
  term00 <- n * (n - 1) * (f11 + f21) * (f12 + f22)
  
  return(term0 + term00)
}

# Simulationsparameter
h1_values <- seq(-1, 1, length.out = 100)  # Werte von h1
h2_values <- seq(-1, 1, length.out = 100)  # Werte von h2
b <- 0.1  # Bandbreite
n <- 10   # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b, n = n))

# Überprüfen, ob NA-Werte vorhanden sind
sum(is.na(grid$z))  # Prüfen auf NA-Werte

# Wenn keine NA-Werte vorhanden sind, Visualisieren der Ergebnisse
Plot1 <- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.1, n = 10)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()



b1 <- 0.05                              # Bandbreite
n1 <- 10                               # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b1, n = n1))

# Visualisieren der Ergebnisse
Plot2<- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.05, n = 10)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()


b0 <- 0.01                               # Bandbreite
n0 <- 10                                # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b0, n = n0))

# Visualisieren der Ergebnisse
Plot0<- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.01, n = 10)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()



b2 <- 0.1                               # Bandbreite
n2 <- 100                               # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b2, n = n2))

# Visualisieren der Ergebnisse
Plot3<- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.1, n = 100)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()


b3 <- 0.05                               # Bandbreite
n3 <- 100                               # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b3, n = n3))

# Visualisieren der Ergebnisse
Plot4<- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.05, n = 100)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()


b00 <- 0.01                               # Bandbreite
n00 <- 100                               # Skalierungsfaktor

# Erstellen eines Datenrasters
grid <- expand.grid(h1 = h1_values, h2 = h2_values)

# Berechnen der Funktion für alle Punkte im Raster
grid$z <- mapply(calculate_expression, grid$h1, grid$h2, 
                 MoreArgs = list(b = b00, n = n00))

# Visualisieren der Ergebnisse
Plot00<- ggplot(grid, aes(x = h1, y = h2, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Visualisierung der Funktion (b = 0.01, n = 100)",
       x = expression(h[1]),
       y = expression(h[2]),
       fill = "Wert") +
  theme_minimal()


library(ggpubr)
combined_plot <- ggarrange(
  Plot1, Plot2, Plot0, Plot3, Plot4, Plot00,
  labels = c("A", "B", "C", "D","E", "F"),  # Labels für die Plots
  ncol = 3, nrow = 2,             # 2x2 Layout
  common.legend = FALSE           # Gemeinsame Legende deaktivieren
)

# Ausgabe anzeigen
print(combined_plot)

################################################################################
################# Erwartungswert Zähler ########################################
################################################################################

library(ggplot2)

# Parameter definieren

phi1 <- 1/3
phi2 <- 1/3
b <- 0.1
sigma <- 1
n <- 100


g <- function(h, b, phi) {
  # Indikatorfunktionen
  ind1 <- (-b / 2 <= h & h <= b / 2)
  ind2 <- (b / 2 - 1 <= h)
  ind3 <- (-1 + b / 2 <= h & h <= -b / 2)
  ind4 <- (h >= -b / 2 & b >= -h + b / 2 & 1 >= -h + b / 2)
  ind5 <- (h <= -b / 2)
  ind6 <- (h >= b / 2)
  ind7 <- (h >= b / 2 & h <= 1 - b / 2)
  
  # Funktionsteil 1
  term1 <- phi * (1 - exp((h - b / 2) / phi)) * (1 - b / 2 + h) * ind1 * ind2
  
  # Funktionsteil 2
  term2 <- phi * (exp((h + b / 2) / phi) - exp((h - b / 2) / phi)) * (1 + h - b / 2) * ind3
  
  # Funktionsteil 3
  term3 <- (phi * (b / 2 - h) - phi^2 * (1 - exp(-(b / 2 - h) / phi))) * ind4
  
  # Funktionsteil 4
  term4 <- phi * exp((h + b / 2) / phi) * 
    (min(1, b / 2 - h) - max(0, -b / 2 - h)) * ind5 +
    phi^2 * (exp(-min(1, b / 2 - h) / phi) - exp(-max(0, -b / 2 - h) / phi)) * ind5
  
  # Funktionsteil 5
  term5 <- phi * exp((-h + b / 2) / phi) *
    (min(1, 1 - h + b / 2) - max(0, 1 - h - b / 2)) * ind6 -
    phi * exp(-1 / phi) * 
    (phi * exp(min(1, 1 - h + b / 2) / phi) - phi * exp(max(0, 1 - h - b / 2) / phi)) * ind6 +
    phi * (exp((-h + b / 2) / phi) - exp((-h - b / 2) / phi)) * (1 - h - b / 2) * ind7
  
  # Funktionsteil 6
  term6 <- (phi * (1 - h - b / 2) + phi^2 * exp(-1 / phi) *
              (1 - exp((1 - h - b / 2) / phi))) * ind1 +
    phi * (1 - exp((-h - b / 2) / phi)) * (1 - h - b / 2) * ind1
  
  # Ergebnis
  result <- term1 + term2 + term3 + term4 + term5 + term6
  return(result)
}


# Gitter von h1 und h2 erstellen

h1_vals <- seq(-1, 1, length.out = 100)
h2_vals <- seq(-1, 1, length.out = 100)

grid <- expand.grid(h1 = h1_vals, h2 = h2_vals)

# Werte berechnen

grid$g_h1 <- g(grid$h1, b, phi1)
grid$g_h2 <- g(grid$h2, b, phi2)

grid$E_Z <- n^2 * sigma^2 * grid$g_h1 * grid$g_h2

# Plot mit ggplot2 erstellen

ggplot(grid, aes(x = h1, y = h2, fill = E_Z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "Visualisierung des Erwartungswert des Zählers",
    x = expression(h[1]),
    y = expression(h[2]),
    fill = expression(Wert)
  ) +
  theme_minimal()


################################################################################
### Plot der wahren Kovarianz###################################################
###############################################################################

library(ggplot2)
library(reshape2)

# Funktion zur Berechnung der exponentiellen Kovarianzmatrix
cov_exponential = function(grid, sigma, phi, method = "difference") {
  if (!method %in% c("difference", "euclidean")) {
    stop("Specify either 'euclidean' or 'difference' method")
  }
  
  if (method == "difference") {
    # Absolute Differenzen berechnen
    dist_matrix_x = abs(outer(grid$x, grid$x, "-"))
    dist_matrix_y = abs(outer(grid$y, grid$y, "-"))
    
    # Kovarianzmatrix berechnen
    covariance = sigma^2 * exp(-dist_matrix_x / phi) * exp(-dist_matrix_y / phi)
  }
  
  return(covariance)
}

# Beispielgitter erzeugen
set.seed(123)
grid = data.frame(x = runif(xdim, 0, 1), y = runif(ydim, 0, 1))

# Parameter für die Kovarianzfunktion
sigma = 1
phi = 3

# Kovarianzmatrix berechnen
cov_matrix = cov_exponential(grid, sigma, phi, method = "difference")

# Matrix in langes Format umwandeln für ggplot
cov_df = melt(cov_matrix)

# Heatmap plotten
ggplot(cov_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = sigma^2 / 2) +
  labs(title = "Exponentielle Kovarianzmatrix (Differenz-Methode)", x = "Index", y = "Index", fill = "Kovarianz") +
  theme_minimal()


