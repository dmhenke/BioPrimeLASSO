#' Gene specific association scores
#'
#' @param gene Gene symbol. e.g. "EGFR"
#' @param network Data frame of parameter-pair scores (range: [0-1])
#'
#' @return A numerical vector named
#' @export
#'
#' @examples
#'
#' head(network)
#'          protein1        protein2 combined_score gene1    gene2
#' 1 ENSP00000000233 ENSP00000356607            0.173  ARF5  RALGPS2
#' 2 ENSP00000000233 ENSP00000427567            0.154  ARF5    FHDC1
#' 3 ENSP00000000233 ENSP00000253413            0.151  ARF5 ATP6V1E1
#' ...
#'
#' scores <- get_scores(gene="EGFR", network)
#' head(scores)
#' ENPEP  LOXL3   SYT1 FCGR2B   AKT1   RARB
#' 0.265  0.197  0.228  0.275  0.792  0.420
get_scores <- function(gene, network){

  tmp <- network[(network$gene1 %in% c(gene) |
                    network$gene2 %in% c(gene)), ]
  tmp <- tmp[tmp$gene1 != "" & tmp$gene2 != "", ]
  tmp <- na.omit(tmp)
  scores <- as.numeric(tmp[["combined_score"]][tmp$gene1 == gene])
  names(scores) <- as.character(tmp[["gene2"]][tmp$gene1 == gene])
  scores[gene] <- 1
  return(scores)
}
