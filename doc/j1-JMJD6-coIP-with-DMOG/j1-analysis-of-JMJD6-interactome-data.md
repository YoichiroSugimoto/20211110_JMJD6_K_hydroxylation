---
title: "j1 Analysis of JMJD6 co-IP data with DMOG treatment"
author: "Yoichiro Sugimoto"
date: "16 December, 2021"
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
   html_document:
     highlight: haddock
     toc: yes
     toc_depth: 2
     keep_md: yes
     fig_width: 5
     fig_height: 5
---


# Environment setup


```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

library("limma")
library("matrixStats")
library("Biostrings")
library("UniProt.ws")

hsdb <- UniProt.ws(9606)

raw.abundance.dt <- fread(
    "../../data/j1-JMJD6-coIP/20210408-JMJD6-interactomics-data.csv"
)

results.dir <- file.path("../../results")
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")

create.dirs(c(
))
```


# Data preprocessing

## Data files


```r
original.sample.names <- colnames(raw.abundance.dt)[3:ncol(raw.abundance.dt)]
sample.names <- original.sample.names %>%
    {gsub("( rep_#|_rep_#)", "_", .)}

setnames(raw.abundance.dt, old = original.sample.names, new = sample.names)

raw.abundance.dt[
  , uniprot_id := str_split_fixed(get("PG.ProteinGroups"), ";", n = 2)[, 1]
]

## Ignore entries without uniprot ID
abundance.dt <- raw.abundance.dt[
  , c("uniprot_id", sample.names), with = FALSE
][!grepl("^ENSEMBL", uniprot_id)]

log2.abundance.dt <- copy(abundance.dt)[
   , (sample.names) := lapply(.SD, log2), .SDcols = sample.names
]

log2.abundance.mat <- as.matrix(log2.abundance.dt[, sample.names, with = FALSE])
rownames(log2.abundance.mat) <- log2.abundance.dt[, uniprot_id]
```

## Sample data

The following samples will be analysed.


```r
sample.dt <- data.table(
    sample_name = sample.names,
    IP_by = str_split_fixed(sample.names, "_", n = 3)[, 1],
    treatment = str_split_fixed(sample.names, "_", n = 3)[, 2],
    replicate = str_split_fixed(sample.names, "_", n = 3)[, 3]
)

print(sample.dt)
```

```
##     sample_name IP_by treatment replicate
##  1:     EV_NT_1    EV        NT         1
##  2:     EV_NT_2    EV        NT         2
##  3:     EV_NT_3    EV        NT         3
##  4:   EV_DMOG_1    EV      DMOG         1
##  5:   EV_DMOG_2    EV      DMOG         2
##  6:   EV_DMOG_3    EV      DMOG         3
##  7:   WTJ6_NT_1  WTJ6        NT         1
##  8:   WTJ6_NT_2  WTJ6        NT         2
##  9:   WTJ6_NT_3  WTJ6        NT         3
## 10: WTJ6_DMOG_1  WTJ6      DMOG         1
## 11: WTJ6_DMOG_2  WTJ6      DMOG         2
## 12: WTJ6_DMOG_3  WTJ6      DMOG         3
## 13:   EDJ6_NT_1  EDJ6        NT         1
## 14:   EDJ6_NT_2  EDJ6        NT         2
## 15:   EDJ6_NT_3  EDJ6        NT         3
## 16: EDJ6_DMOG_1  EDJ6      DMOG         1
## 17: EDJ6_DMOG_2  EDJ6      DMOG         2
## 18: EDJ6_DMOG_3  EDJ6      DMOG         3
```

## Extraction of max K ratio

Maximum K ratio of proteins by 10-mer sliding window will be extracted.
Note that genes with ENSEMBL:xx ids were ignored. 


```r
all.protein.feature.per.pos.dt <- file.path(
    j0.res.dir, "all_protein_feature_per_position.csv"
) %>% fread

maxK.per.protein.dt <- all.protein.feature.per.pos.dt[
  , list(max_k_ratio = max(K_ratio), aa_len = max(position)),
    by = list(Accession, uniprot_id)
]

maxK.per.protein.dt <- maxK.per.protein.dt[aa_len > 20]
```

## Extraction of gene name



```r
gene.name.dt <- log2.abundance.dt[, uniprot_id] %>%
    {select(
         hsdb,
         keys = .,
         keytype = "UNIPROTKB",
         columns = c("GENES"),
         multiVals = "list"
     )} %>%
    data.table 
```

```
## Uniprot limits queries with a large amount of keys. It's recommended that the select method be invoked with fewer than 100 keys or the query may fail.
```

```
## Getting extra data for A0A075B6S2, A8MWD9, E9PAV3... (400 total)
```

```
## Getting extra data for Q2TAY7, Q2TBE0, Q3B726... (173 total)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```r
gene.name.dt[
  , gene_name := str_split_fixed(
        GENES, " ", n = 2
    )[, 1]]

setnames(gene.name.dt, old = "UNIPROTKB", new = "uniprot_id")

maxK.per.protein.dt <- merge(
    gene.name.dt[, .(uniprot_id, gene_name)],
    maxK.per.protein.dt,
    by = "uniprot_id",
    all.x = TRUE
) %>%
    {.[!is.na(Accession)]}
```


# Heatmap analysis

Proteins with different abundance as a function of JMJD6 IP or DMOG will be identified, and then these proteins will be analysed by the heatmap.

In detail, the following comparisons will be performed:

- DMOGvsNAbyJMJD6: IPed by JMJD6 with DMOG or without treatment
- JMJD6vsEDJMJD6inNA: JMJD6 or EDJMJD6 IP without treatment
- JMJD6vsEDJMJD6inDMOG: JMJD6 or EDJMJD6 IP with DMOG treatment
- JMJD6vsNAinNA: JMJD6 or no antibody IP without treatment
- JMJD6vsNAinDMOG: JMJD6 or no antibody IP with DMOG treatment

Statistical test for all these comparisons will be performed, and all genes identified as significant (FDR < 0.01) by any comparison will be analysed by heatmap.

## Statistical test


```r
exp.condition <-
    colnames(log2.abundance.mat) %>%
    {paste(
        str_split_fixed(., "_", n = 3)[, 1],
        str_split_fixed(., "_", n = 3)[, 2],
        sep = "."
     )} %>%
    factor(levels = c(
               "WTJ6.NT", "WTJ6.DMOG",
               "EDJ6.NT", "EDJ6.DMOG",
               "EV.NT", "EV.DMOG"
           ))

design <- model.matrix(~ 0 + exp.condition)
colnames(design) <- levels(exp.condition)

fit <- lmFit(log2.abundance.mat, design)

contrast.mat <- makeContrasts(
    DMOGvsNAbyJMJD6 = WTJ6.DMOG - WTJ6.NT,
    JMJD6vsEDJMJD6inNA = WTJ6.NT - EDJ6.NT,
    JMJD6vsEDJMJD6inDMOG = WTJ6.DMOG - EDJ6.DMOG,
    JMJD6vsNAinNA = WTJ6.NT - EV.NT,
    JMJD6vsNAinDMOG = WTJ6.DMOG - EV.DMOG,    
    levels = design
)

fit2 <- contrasts.fit(fit, contrast.mat)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- decideTests(fit2, method = "global")

exportTopTable <- function(coef, fit2){
    res.dt <- data.table(
        topTable(fit2, coef = coef, number = Inf, sort = "none"),
        keep.rownames = TRUE
    )

    res.dt <- res.dt[, c("rn", "AveExpr", "logFC", "adj.P.Val"), with = FALSE]
    
    setnames(
        res.dt,
        old = colnames(res.dt),
        new = c("uniprot_id", paste0(coef, "_", c("AveExpr", "logFC", "adj.P.Val")))
    )

    setkey(res.dt)

    return(res.dt)
}


res.dts <- lapply(
    colnames(contrast.mat),
    exportTopTable,
    fit2 = fit2
)

res.dt <- Reduce(function(...) merge(..., all = TRUE), res.dts)

fdrs <- do.call(pmin, res.dt[, grep("_adj.P.Val", colnames(res.dt)), with = FALSE])
sig.genes <- res.dt[fdrs < 0.01, uniprot_id]

sig.protein.mat <- log2.abundance.mat[rownames(log2.abundance.mat) %in% sig.genes, ]

col.orders <- sample.dt[
    order(match(treatment, c("NT", "DMOG")))
][
    order(match(IP_by, c("EV", "WTJ6", "EDJ6")))
][, sample_name]

sig.protein.mat <- sig.protein.mat[, col.orders]

rownames(sig.protein.mat) %>%
    {.[!(. %in% maxK.per.protein.dt[, uniprot_id])]}
```

```
## character(0)
```

```r
nm.conv <- setNames(
    maxK.per.protein.dt[, gene_name],
    nm = maxK.per.protein.dt[, uniprot_id]
)

rownames(sig.protein.mat) <- nm.conv[rownames(sig.protein.mat)]
```

## Plot heatmap


```r
library("pheatmap")

cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
}

scaled.fc.mat <- apply(sig.protein.mat, 1, cal_z_score)

c <- cor(scaled.fc.mat, method = "spearman") 
d <- as.dist(1-c)

pheatmap(
    t(scaled.fc.mat),
    cluster_cols = FALSE,
    clustering_distance_rows = d,
    clustering_method = "ward.D2",
    border_color = FALSE,
    show_rownames = TRUE
)
```

![](j1-analysis-of-JMJD6-interactome-data_files/figure-html/heatmap2-1.png)<!-- -->

# DMOG dependent interactome



```r
res.summry.dt <- merge(
    maxK.per.protein.dt,
    res.dt,
    by = "uniprot_id",
    all.y = TRUE
)

fwrite(
    res.summry.dt,
    file.path(results.dir, "differential-protein-abundance.csv")
)
```


# K ratio



```r
m.res.summry.dt <- melt(
    res.summry.dt,
    measure.vars = grep("_logFC$", colnames(res.summry.dt), value = TRUE),
    id.vars = c("uniprot_id", "Accession", "max_k_ratio"),
    value.name = "log2FC",
    variable.name = "comparison_name"
) %>%
    {.[, `:=`(
         IP_name = str_split_fixed(comparison_name, "_", n = 2)[, 1])
       ]} 

m.res.summry.dt[, IP_name := factor(
               IP_name, levels = c(
                            "DMOGvsNAbyJMJD6", "JMJD6vsNAinNA", "JMJD6vsNAinDMOG",
                            "JMJD6vsEDJMJD6inNA", "JMJD6vsEDJMJD6inDMOG"
                        )
           )][, `:=`(
                   treated_with = str_split_fixed(
                       IP_name, "in", n = 2
                   )[, 2] %>%
                       {gsub("NA", "no treatment", .)} %>%
                       factor(levels = c("no treatment", "DMOG")),
                   max_k_ratio_md = case_when(
                       max_k_ratio <= 0.2 ~ "0.1-0.2",
                       TRUE ~ as.character(max_k_ratio)
                   )
               )]

print("Protein number by k ratio:")
```

```
## [1] "Protein number by k ratio:"
```

```r
m.res.summry.dt[
    comparison_name == "DMOGvsNAbyJMJD6_logFC",
    table(max_k_ratio)
]
```

```
## max_k_ratio
## 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9   1 
##   9  60 199 174  59  33  18  14   7   1
```

```r
sig.th <- 0.05

sig.dt <- m.res.summry.dt[comparison_name == "DMOGvsNAbyJMJD6_logFC"] %$%
    pairwise.wilcox.test(
        x = log2FC,
        g = max_k_ratio_md,
        p.adjust.method = "none"
    )$p.value[, 1] %>%
    stack %>% data.table %>%
    {.[, padj := p.adjust(values, method = "holm")]} %>%
    {.[, `:=`(
         sig_mark = case_when(
             padj < sig.th * 0.1 ~ "**",
             padj < sig.th ~ "*",
             TRUE ~ NA_character_
         ),
         max_k_ratio_md = ind
     )]}

merge(
    m.res.summry.dt[comparison_name == "DMOGvsNAbyJMJD6_logFC"],
    sig.dt,
    by = "max_k_ratio_md", all.x = TRUE
) %>%
    ggplot(
        aes(
            x = max_k_ratio_md,
            y = log2FC,
            color = max_k_ratio != 1
        )
    ) +
    geom_hline(yintercept = 0, color = "gray60") +
    ## ggbeeswarm::geom_quasirandom(color = "gray60") +
    geom_boxplot(notch = FALSE, fill = "white", alpha = 0.6, outlier.shape = NA) +
    stat_summary(
        geom = 'text', aes(label = sig_mark),
        fun = function(x){boxplot.stats(x)$stats[5]}, 
        vjust = -0.8, color = "black", size = 5
    ) +
    coord_cartesian(ylim = c(-5, 5)) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray80")) +
    ## viridis::scale_fill_viridis(discrete = TRUE, direction = -1) +
    xlab("Maximum K ratio in 10 mer sliding window") +
    ylab("Protein abundance log2 FC by JMJD6 IP\n(DMOG treatment vs no treamtent)") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        aspect.ratio = 0.8
    )
```

```
## Warning: Removed 7 rows containing missing values (geom_text).
```

![](j1-analysis-of-JMJD6-interactome-data_files/figure-html/k ratio-1.png)<!-- -->

# Session information



```r
sessioninfo::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.0.0 (2020-04-24)
##  os       CentOS Linux 7 (Core)       
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_GB.UTF-8                 
##  ctype    en_GB.UTF-8                 
##  tz       Europe/London               
##  date     2021-12-16                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package       * version  date       lib source        
##  AnnotationDbi   1.50.0   2020-04-27 [1] Bioconductor  
##  assertthat      0.2.1    2019-03-21 [1] CRAN (R 4.0.0)
##  Biobase         2.48.0   2020-04-27 [1] Bioconductor  
##  BiocFileCache   1.12.0   2020-04-27 [1] Bioconductor  
##  BiocGenerics  * 0.34.0   2020-04-27 [1] Bioconductor  
##  Biostrings    * 2.56.0   2020-04-27 [1] Bioconductor  
##  bit             1.1-15.2 2020-02-10 [1] CRAN (R 4.0.0)
##  bit64           0.9-7    2017-05-08 [1] CRAN (R 4.0.0)
##  bitops          1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
##  blob            1.2.1    2020-01-20 [1] CRAN (R 4.0.0)
##  cachem          1.0.5    2021-05-15 [1] CRAN (R 4.0.0)
##  cli             3.0.1    2021-07-17 [1] CRAN (R 4.0.0)
##  colorspace      1.4-1    2019-03-18 [1] CRAN (R 4.0.0)
##  crayon          1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
##  curl            4.3      2019-12-02 [1] CRAN (R 4.0.0)
##  data.table    * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)
##  DBI             1.1.0    2019-12-15 [1] CRAN (R 4.0.0)
##  dbplyr          1.4.4    2020-05-27 [1] CRAN (R 4.0.0)
##  digest          0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
##  dplyr         * 1.0.0    2020-05-29 [1] CRAN (R 4.0.0)
##  ellipsis        0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
##  evaluate        0.14     2019-05-28 [1] CRAN (R 4.0.0)
##  farver          2.0.3    2020-01-16 [1] CRAN (R 4.0.0)
##  fastmap         1.0.1    2019-10-08 [1] CRAN (R 4.0.0)
##  generics        0.0.2    2018-11-29 [1] CRAN (R 4.0.0)
##  ggplot2       * 3.3.1    2020-05-28 [1] CRAN (R 4.0.0)
##  glue            1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
##  gtable          0.3.0    2019-03-25 [1] CRAN (R 4.0.0)
##  htmltools       0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
##  httr            1.4.1    2019-08-05 [1] CRAN (R 4.0.0)
##  IRanges       * 2.22.1   2020-04-28 [1] Bioconductor  
##  khroma        * 1.3.0    2019-10-26 [1] CRAN (R 4.0.0)
##  knitr         * 1.28     2020-02-06 [1] CRAN (R 4.0.0)
##  labeling        0.3      2014-08-23 [1] CRAN (R 4.0.0)
##  lifecycle       0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
##  limma         * 3.44.1   2020-04-28 [1] Bioconductor  
##  magrittr      * 1.5      2014-11-22 [1] CRAN (R 4.0.0)
##  matrixStats   * 0.56.0   2020-03-13 [1] CRAN (R 4.0.0)
##  memoise         2.0.0    2021-01-26 [1] CRAN (R 4.0.0)
##  munsell         0.5.0    2018-06-12 [1] CRAN (R 4.0.0)
##  pheatmap      * 1.0.12   2019-01-04 [1] CRAN (R 4.0.0)
##  pillar          1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
##  pkgconfig       2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
##  purrr           0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
##  R6              2.4.1    2019-11-12 [1] CRAN (R 4.0.0)
##  rappdirs        0.3.1    2016-03-28 [1] CRAN (R 4.0.0)
##  RColorBrewer    1.1-2    2014-12-07 [1] CRAN (R 4.0.0)
##  Rcpp            1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
##  RCurl         * 1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
##  rlang           0.4.11   2021-04-30 [1] CRAN (R 4.0.0)
##  rmarkdown     * 2.2      2020-05-31 [1] CRAN (R 4.0.0)
##  RSQLite       * 2.2.0    2020-01-07 [1] CRAN (R 4.0.0)
##  S4Vectors     * 0.26.0   2020-04-27 [1] Bioconductor  
##  scales          1.1.1    2020-05-11 [1] CRAN (R 4.0.0)
##  sessioninfo     1.1.1    2018-11-05 [1] CRAN (R 4.0.5)
##  statmod         1.4.34   2020-02-17 [1] CRAN (R 4.0.0)
##  stringi         1.4.6    2020-02-17 [1] CRAN (R 4.0.0)
##  stringr       * 1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
##  tibble          3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
##  tidyselect      1.1.0    2020-05-11 [1] CRAN (R 4.0.0)
##  UniProt.ws    * 2.28.0   2020-04-27 [1] Bioconductor  
##  vctrs           0.3.1    2020-06-05 [1] CRAN (R 4.0.0)
##  withr           2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
##  xfun            0.14     2020-05-20 [1] CRAN (R 4.0.0)
##  XVector       * 0.28.0   2020-04-27 [1] Bioconductor  
##  yaml            2.2.1    2020-02-01 [1] CRAN (R 4.0.0)
##  zlibbioc        1.34.0   2020-04-27 [1] Bioconductor  
## 
## [1] /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/hydroxylation_by_JMJD6/lib/R/library
```
