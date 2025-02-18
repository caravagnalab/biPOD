#' Xenografts Mouse Model Dataset
#'
#' A dataset containing tumor volume measurements from xenograft mouse models under
#' different treatment conditions. The data tracks tumor growth in treated and
#' untreated mice over time, with measurements beginning at treatment initiation.
#'
#' @format A data frame with 159 rows and 4 columns:
#' \describe{
#'   \item{mouse}{Unique identifier for each mouse in the study}
#'   \item{tumour_volume}{Measured volume of the tumor}
#'   \item{time}{Time in days since treatment initiation (day 0 represents treatment start)}
#'   \item{treatment_group}{Binary variable indicating treatment status (treated/untreated)}
#' }
#'
#' @source Data from Sauer et al. https://doi.org/10.15252/emmm.202215729
#'
"xenografts"
