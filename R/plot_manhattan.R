#' Manhattan plots of bplasso results
#'
#' @param gene Gene symbol. e.g. "EGFR"
#' @param resIn Path to "results_omic".RData file.
#' @param subplotChr Chromosome number for subplot
#' @param dependency Dependency data frame which includes gene column
#' @param gene_info Gene information data frame. Columns should include : hgnc_symbol, chromosome_name, start_position
#' @param dir_save directory in which plots will be saved
#'
#' @return Three Manhattan plots are written into the dir_save folder
#' @export
#'
#' @examples
#' Gene <- "EGFR"
#' results_EGFR <- "../Outputs/EGFR_demeter2_CNV.RData"
#' eleven <- 11
#' D2 <- demeter2
#' dir_save <- "../Outputs/Graphics/"
#' plot_manhattan(gene=Gene,resIn=results_EGFR,subplotChr=eleven,dependency=D2,dir_save=fldr_plts)


plot_manhattan <- function(gene,resIn,subplotChr=NA,dependency,gene_info,dir_save){
  load(resIn)
  y <- dependency[, gene]
  y <- y[!is.na(y)]
  correl <- results_omic$cor2score
  aframe <- data.frame(
    gene = names(correl),
    gene_info[match(names(correl), gene_info$hgnc_symbol), ],
    correl
  )
  aframe <- aframe[order(aframe$chromosome_name, aframe$start_position),]
  aframe$rank <- 1:nrow(aframe)

  subm <- results_omic
  subm_betas <- subm$betas[match(aframe$gene, rownames(subm$betas)),]
  subm_betas$betas_pen[which(subm_betas$betas_pen==0)]<- NA;subm_betas$betas[which(subm_betas$betas==0)]<- NA;
  aframe$betas_pen <- subm_betas$betas_pen
  aframe$betas <- subm_betas$betas
  aframe$betalogic <- apply(cbind(aframe$betas_pen,aframe$betas),1,function(x){
    if(is.na(x[1])&is.na(x[2])) NA else if(!is.na(x[1])&is.na(x[2])) 'beta_pen' else if(!is.na(x[1])&!is.na(x[2])) 'both' else 'beta'
  })
  aframe <- aframe[!is.na(aframe$chromosome_name), ]
  aframe$betasize <- aframe$betas;  aframe$betasize[!is.na(aframe$betas_pen)] <- aframe$betas_pen[!is.na(aframe$betas_pen)]
  aframe$betasize <- abs(aframe$betasize)
  # labs
  lab_x <- "Gene ordered by genomic coordinate"
  lab_y <- paste0("Correlation (r)\n",gene," Dependency with 'omic")

  g_full <- ggplot(aframe,
                   aes(rank, correl, color = chromosome_name)) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(size=0.3) +
    scale_color_manual(guide='none', breaks=c(1:22,"beta_pen","both","beta"),
                       values = c(rep(c("black", "grey"), 11),"darkblue","purple","darkred")) +
    labs(x=lab_x,y=lab_y)+
    theme_classic()

  g_fullLab <- g_full +
    ggrepel::geom_text_repel(
      min.segment.length = 0,
      force=1,direction='both',max.overlaps=100,
      max.time = .3, max.iter = 1e5,
      data =aframe[!is.na(aframe$betalogic),],
      aes(label = gene, color = betalogic,size=betasize))+
    scale_size(guide = 'none')
  ggsave(filename = paste0(dir_save,"/PlotManhattan_",gene,".pdf"),plot = g_full,width = 10,height = 5)
  ggsave(filename = paste0(dir_save,"/PlotManhattan_",gene,"_labs.pdf"),plot = g_fullLab,width = 10,height = 5)

  # Subplot of chr
  if(is.na(subplotChr)){
    which_chrome <- aframe[which(aframe$gene==gene),"chromosome_name"]
  } else which_chrome <- subplotChr
  which_chromePos <- which(aframe$chromosome_name==which_chrome)
  # x min max
  which_chrome_minX <- min(which_chromePos)
  which_chrome_maxX <- max(which_chromePos)
  # y min max
  which_chrome_minY <-  min (aframe[which(aframe$chromosome_name==which_chrome),'correl'])
  which_chrome_maxY <-  max (aframe[which(aframe$chromosome_name==which_chrome),'correl'])

  g_chrGene <- g_fullLab+
    coord_cartesian(xlim=c(which_chrome_minX,which_chrome_maxX),
                    ylim=c(which_chrome_minY,which_chrome_maxY))+
    expand_limits(x = 0, y = 0)+
    scale_x_continuous(expand = c(0, 0))+ scale_y_continuous(expand = c(0, 0))+
    labs(x=paste0("Chromosome ",which_chrome))
  ggsave(filename = paste0(dir_save,"/PlotManhattan_",gene,"_chr",which_chrome,".pdf"),plot = g_chrGene,width = 4,height = 2)
  plot(g_fullLab)
}
