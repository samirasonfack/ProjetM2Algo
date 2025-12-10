#' @title Compare deux séries temporelles avant et après alignement DTW (version simplifiée)
#'
#' @description
#' Applique l'algorithme Dynamic Time Warping (DTW) avec la contrainte de la
#' bande de Sakoe-Chiba pour aligner deux séries temporelles numériques X et Y,
#' et affiche les séries brutes et les séries alignées côte à côte.
#'
#' @param X La première série temporelle (vecteur numérique).
#' @param Y La seconde série temporelle (vecteur numérique).
#' @param radius Demi-largeur de la bande de Sakoe-Chiba (contrainte DTW).
#' @param dtw_func La fonction Rcpp du DTW avec bande à utiliser (ex: dtw_sakoe_chiba_rcpp).
#'
#' @return Un objet de type 'grob' (grid graphical object) contenant les deux graphiques
#'         organisés côte à côte.
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal scale_color_manual
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate # Importation spécifique de mutate pour dplyr
#' @importFrom gridExtra grid.arrange
#' @export
plot_dtw_comparison_simple <- function(X, Y, radius, dtw_func) {

  # --- 1. Préparation pour le graphique 'Avant DTW' (Séries Brutes) ---

  # Créer des dataframes pour les séries brutes
  df_X <- data.frame(Index = seq_along(X), Value = X, Series = "Série X")
  df_Y <- data.frame(Index = seq_along(Y), Value = Y, Series = "Série Y")
  df_before <- rbind(df_X, df_Y)

  plot_before <- ggplot2::ggplot(df_before, ggplot2::aes(x = Index, y = Value, color = Series)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Avant DTW : Séries Brutes",
                  y = "Valeur",
                  color = "Série") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("Série X" = "blue", "Série Y" = "red"))

  # --- 2. Calcul et Alignement DTW ---
  message(paste("Calcul DTW avec radius =", radius, "..."))

  # Utilisation de la fonction Rcpp passée en argument
  dtw_result <- dtw_func(X, Y, radius = radius)

  # Extraire le chemin d'alignement (indices)
  path <- dtw_result$path
  idx_x <- sapply(path, `[[`, 1)
  idx_y <- sapply(path, `[[`, 2)

  aligned_X <- X[idx_x]
  aligned_Y <- Y[idx_y]

  df_aligned <- data.frame(
    Index = seq_along(aligned_X),
    X_Aligned = aligned_X,
    Y_Aligned = aligned_Y
  )

  # Réorganiser pour ggplot (format long) en utilisant tidyr::pivot_longer
  df_after_long <- df_aligned %>%
    tidyr::pivot_longer(
      cols = c(X_Aligned, Y_Aligned),
      names_to = "Series",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Series = ifelse(Series == "X_Aligned", "Série X", "Série Y"))


  # --- 3. Génération du graphique 'Après DTW' (Séries Alignées) ---

  plot_after <- ggplot2::ggplot(df_after_long, ggplot2::aes(x = Index, y = Value, color = Series)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Après DTW (Dist. =", round(dtw_result$distance, 2), ", Radius =", radius, ")"),
                  x = paste("Indice Alignement (Longueur :", length(aligned_X), ")"),
                  y = "Valeur (Alignée)",
                  color = "Série") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("Série X" = "blue", "Série Y" = "red"))

  # --- 4. Retour des graphiques organisés ---

  gridExtra::grid.arrange(plot_before, plot_after, ncol = 2)
}
