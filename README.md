## Bio-primed machine learning to enhance discovery of relevant biomarkers

<p align="middle">
  <img src="BioPrimeLASSO_overview.png" width="75%"/>
</p>

### Introduction

Precision medicine relies on identifying reliable biomarkers for gene dependencies to tailor individualized therapeutic strategies. The advent of high-throughput technologies presents unprecedented opportunities to explore molecular disease mechanisms but also challenges due to high dimensionality and collinearity among features. Traditional statistical methods often fall short in this context, necessitating novel computational approaches that harness the full potential of big data in bioinformatics. Here, we introduce a novel machine learning approach extending the Least Absolute Shrinkage and Selection Operator (LASSO) regression framework to incorporate biological knowledge, such as protein-protein interaction databases, into the regularization process. This bio-primed approach prioritizes variables that are both statistically significant and biologically relevant. Applying our method to multiple dependency datasets, we identified biomarkers which traditional methods overlooked. Our biologically informed LASSO method effectively identifies relevant biomarkers from high-dimensional collinear data, bridging the gap between statistical rigor and biological insight. This method holds promise for advancing personalized medicine by uncovering novel therapeutic targets and understanding the complex interplay of genetic and molecular factors in disease.

------------------------------------------
### Reproducibility
Analysis code to reproduce results described in our manuscript can be found [here](https://github.com/dmhenke/BioPrimeLASSO/tree/main/Reproducibility).

------------------------------------------

### R Package Walkthrough

#### 1)  Installation

Our R package called **BioPrimeLASSO** requires the following R packages to be installed: [glmnet](https://github.com/cran/glmnet) and [ggplot2](https://github.com/tidyverse/ggplot2).

``` r
install.packages("devtools")
devtools::install_github("dmhenke/BioPrimeLASSO")
```

#### 2)  Load toy data (total size ~20Mb)

  In this toy example we will use BioPrimeLASSO to discover copy number biomarkers for _EGFR_ dependency. BioPrimeLASSO also makes use of Protein-Protein interaction information from STRING DB. Please download the following three files:

1. Copy number variation ([cnv_EGFR.tsv](https://drive.google.com/file/d/1aqWQcxg3CgFGSCrpcz6ElZ1NX1t6u27U/view?usp=drive_link))
2. Dependency data ([demeter2_EGFR.tsv](https://drive.google.com/file/d/13gyAJg6XHofzWbMuNtSEECR69WLPpvwG/view?usp=drive_link))
3. Protein-protein interaction network ([ppi_w_symbols_EGFR.tsv](https://drive.google.com/file/d/1npIekQYq_GgpyF6z2NLUMgIoLtUaaQ9M/view?usp=drive_link))

``` r
cnv <- read.csv("./cnv_EGFR.tsv",sep = '\t',header=T)
ppi <- read.csv("./ppi_w_symbols_EGFR.tsv",sep = '\t',header=T)
demeter2 <- read.csv("./demeter2_EGFR.tsv",sep = '\t',header=T)
```

#### 2.1) Load supplemental information (optional)
Next, we load some information for each gene including genomic location using the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) R package.
``` r
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = colnames(cnv),
  mart = mart)

chrs <- as.character(1:22)
gene_info <- gene_info[gene_info$chromosome_name %in% chrs, ]
uniq <- names(which(table(gene_info$hgnc_symbol) == 1))
gene_info <- gene_info[gene_info$hgnc_symbol %in% uniq, ]
gene_info$chromosome_name <- factor(
  gene_info$chromosome_name, levels = chrs)
```

#### 3)  Define gene of interest: EGFR

``` r
GoI <- "EGFR"
```

#### 4)  Setup data objects for analysis

``` r
# Dependency score resource: demeter2
y <- demeter2[,GoI]
names(y) <- rownames(demeter2)

# Identify 'omic information to test against dependency score: cnv
X_omic <- cnv

## Refine population to overlapping cell lines
ok_cells <- intersect(names(y), rownames(X_omic))
X_omic_OK  <- X_omic[ok_cells, ]
y_ok <- y[ok_cells]

## Remove features without variance ####
X_omic_OK <- X_omic_OK[, apply(X_omic_OK, 2, var) > 0]

### Generate scores
# Format: colnames(network) <- c("combined_score","gene1","gene2")
scores <- get_scores(gene=GoI, network=ppi)
```

#### 5)  Run BioPrimeLASSO

``` r
results_omic <- bplasso(
  scale(X_omic_OK), y_ok, scores,
  n_folds = 10, phi_range = seq(0, 1, length = 30))

# Add Pearson correlation: cor2score
results_omic$cor2score <- cor(
  X_omic_OK, y_ok,
  use = "pairwise.complete")[,1]

# Save results
file_results <- paste0("./",GoI,"_demeter2_CNV.RData")
save(results_omic,file = file_results)
```

#### 6)  Visualize results

``` r
## Correlation of Dependency score and CNV for each gene overlaying bio-primed LASSO & baseline LASSO hits
plot_manhattan(gene=GoI,
  resIn=file_results,
  subplotChr=11,
  dependency=demeter2,
  gene_info=gene_info,
  dir_save="./")
```

-----------------------------------------------------------

### Data

For full analysis and to reproduce the results in our manuscript please use the following files (total size ~2Gb):

1. Protein-protein interaction network ([ppi_w_symbols.tsv](https://drive.google.com/file/d/1-Flap0yM1Ba4d8ibVYs6ha82snsmAu-v/view?usp=drive_link))
2. Copy number variation ([cnv.tsv](https://drive.google.com/file/d/1dtKIOnx_lVn5glp67ItjPbiSdE10ZFFm/view?usp=drive_link))
3. RNA expression ([rna.tsv](https://drive.google.com/file/d/1oNrQNUHXkjVy0HgNnNst9yR_cqKFLYyZ/view?usp=drive_link))
4. Demeter2 dependency data ([demeter2.tsv](https://drive.google.com/file/d/1loo9kdMwAUYoJrBCwe3Dk1b9TDDyY72e/view?usp=drive_link))
5. Chronos dependency data ([chronos.tsv](https://drive.google.com/file/d/1_TvvBO26EFDXR7nIXUmEYt919FaD2Ve1/view?usp=drive_link))

These files were originally downloaded from the DepMap [webportal](https://depmap.org/portal/) (22Q2) and STRING DB [website](https://string-db.org/).

-----------------------------------------------------------
## License

`BioPrimeLASSO` uses GNU General Public License GPL-3.

-----------------------------------------------------------

