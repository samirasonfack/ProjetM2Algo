
#' DTW algorithm testing with a case example
x <- c(1, 3, 4, 9,2)
y <- c(1, 3, 7, 8, 9)
dtw_result <- dtw(x, y)

dtw_result



# Vecteurs de tailles à tester
taille_vecteurs <- seq(50, 1000, by = 50)

# Stocker les temps
temps <- numeric(length(taille_vecteurs))

# Boucle de test
for (k in seq_along(taille_vecteurs)) {
  n <- taille_vecteurs[k]
  x <- runif(n)  # série aléatoire
  y <- runif(n)

  temps[k] <- system.time(dtw(x, y))["elapsed"]
}

# Regrouper les résultats
resultats <- data.frame(
  Taille = taille_vecteurs,
  Temps = temps
)

plot(resultats$Taille, resultats$Temps, type = "b", pch = 19, col = "blue",
     main = "Complexité expérimentale de l'algorithme DTW",
     xlab = "Taille des vecteurs (n)",
     ylab = "Temps d'exécution (secondes)")

plot(log(resultats$Taille), log(resultats$Temps), type = "b", pch = 19, col = "blue",
     main = "Complexité expérimentale de l'algorithme DTW",
     xlab = "Taille des vecteurs (n)",
     ylab = "Temps d'exécution (secondes)")

