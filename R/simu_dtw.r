#' Simulation du temps d'exécution de l'algorithme DTW
#'
#' Cette fonction mesure expérimentalement le temps d'exécution
#' de l'algorithme DTW (Dynamic Time Warping) pour des vecteurs de
#' différentes tailles. Elle permet de visualiser la complexité
#' expérimentale de l'algorithme et d'estimer la pente sur un
#' graphique log-log.
#'
#' @param debut Nombre. Taille minimale des vecteurs simulés (par défaut 400).
#' @param fin Nombre. Taille maximale des vecteurs simulés (par défaut 5000).
#' @param pas Nombre. Pas entre chaque taille de vecteur (par défaut 50).
#' @param afficher_plots Booléen. Indique si les graphiques doivent être affichés (par défaut TRUE).
#'
#' @return Une liste contenant :
#' \describe{
#'   \item{resultats}{Un data.frame avec la taille des vecteurs et les temps d'exécution.}
#'   \item{regression}{L'objet lm de la régression log-log.}
#'   \item{pente}{La pente estimée de la régression log-log.}
#'   \item{intercept}{L'ordonnée à l'origine de la régression log-log.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- simulation_dtw(debut = 100, fin = 1000, pas = 100)
#' }
#' @importFrom stats runif na.omit lm coef
#' @importFrom graphics abline
#' @export
simulation_dtw <- function(
    debut = 400,
    fin = 5000,
    pas = 50,
    afficher_plots = TRUE
) {
  # générer les tailles de vecteurs
  taille_vecteurs <- seq(debut, fin, by = pas)

  # vecteur vide pour stocker les temps
  temps <- numeric(length(taille_vecteurs))

  # boucle de calcul des temps d'exécution
  for (k in seq_along(taille_vecteurs)) {
    n <- taille_vecteurs[k]
    x <- runif(n)
    y <- runif(n)

    temps[k] <- system.time(dtw(x, y))["elapsed"]
  }

  # rassembler les résultats
  resultats <- data.frame(
    Taille = taille_vecteurs,
    Temps = temps
  )

  # --- gérer les temps nuls ou négatifs pour éviter log(0) ---
  resultats$Temps[resultats$Temps <= 0] <- NA
  resultats <- na.omit(resultats)

  if (afficher_plots) {
    # --- Graphique taille vs temps ---
    plot(resultats$Taille, resultats$Temps, type = "b", pch = 19, col = "blue",
         main = "Complexité expérimentale de l'algorithme DTW",
         xlab = "Taille des vecteurs (n)",
         ylab = "Temps d'exécution (secondes)")

    # --- Graphique log-log ---
    plot(log(resultats$Taille), log(resultats$Temps), type = "b", pch = 19, col = "blue",
         main = "Complexité expérimentale (log-log) de DTW",
         xlab = "log(Taille des vecteurs (n))",
         ylab = "log(Temps d'exécution (secondes))")
  }

  # --- Régression log-log ---
  regression <- lm(log(Temps) ~ log(Taille), data = resultats)
  pente <- coef(regression)[2]
  intercept <- coef(regression)[1]

  # afficher la pente estimée
  cat("Pente estimée sur le graphique log-log :", round(pente, 2), "\n")

  # ajouter la droite de régression si les plots ont été générés
  if (afficher_plots) {
    abline(regression, col = "red", lwd = 2)
  }

  # retourner les résultats utiles
  return(list(
    resultats = resultats,
    regression = regression,
    pente = pente,
    intercept = intercept
  ))
}




