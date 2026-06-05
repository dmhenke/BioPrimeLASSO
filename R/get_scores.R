#' Gene-specific association scores from a PPI network
#'
#' Extracts and returns the protein-protein interaction (PPI) scores for all
#' network neighbors of a given gene. The query gene itself is assigned a score
#' of 1 (maximum relevance).
#'
#' @param gene A character string giving the gene symbol of interest (e.g., \code{"EGFR"}).
#' @param network A data frame of protein-pair interaction scores with the following
#'   required columns:
#'   \describe{
#'     \item{protein1}{Ensembl protein ID for partner 1.}
#'     \item{protein2}{Ensembl protein ID for partner 2.}
#'     \item{combined_score}{Numeric interaction confidence score, range [0, 1].}
#'     \item{gene1}{HGNC symbol for partner 1.}
#'     \item{gene2}{HGNC symbol for partner 2.}
#'   }
#'   Typically sourced from STRING DB and formatted accordingly.
#'
#' @return A named numeric vector of interaction scores for all neighbors of
#'   \code{gene} in the network, with the query gene set to 1.
#' @export
#'
#' @examples
#' \dontrun{
#'   # network columns: protein1, protein2, combined_score, gene1, gene2
#'   scores <- get_scores(gene = "EGFR", network = ppi)
#'   head(scores)
#'   # ENPEP  LOXL3   SYT1 FCGR2B   AKT1   RARB
#'   # 0.265  0.197  0.228  0.275  0.792  0.420
#' }
get_scores <- function(gene, network){

  tmp <- network[(network$gene1 %in% c(gene) |
                    network$gene2 %in% c(gene)), ]
  tmp <- tmp[tmp$gene1 != "" & tmp$gene2 != "", ]
  tmp <- stats::na.omit(tmp)
  scores <- as.numeric(tmp[["combined_score"]][tmp$gene1 == gene])
  names(scores) <- as.character(tmp[["gene2"]][tmp$gene1 == gene])
  scores[gene] <- 1
  return(scores)
}
