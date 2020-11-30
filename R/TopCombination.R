#' @title Combination Drug Selection.
#' @description  Select drug combinations by combination therapeutic score and FDR of combination therapeutic score.
#' @details Input raw drug combination result and return the top drug combinations.
#' @param Drug.combination raw drug combination result from DrugCombination function.
#' @param Combination.FDR The FDR threshold to select drug combination. The default value is 0.1.
#' @param Min.combination.score The Combination therapeutic score threshold to select drug combination. The default value is 1.
#' @return A data frame of selected drug combinations.
#' @export


TopCombination <- function(Drug.combination=Drug.combinations,
                           Combination.FDR=0.1,
                           Min.combination.score=1
){
  Drug.combination <- subset(Drug.combination, Combination.therapeutic.score > Min.combination.score & FDR < Combination.FDR)
  Drug.combination <- Drug.combination[order(Drug.combination$Combination.therapeutic.score, decreasing = T),]
  return(Drug.combination)
}
