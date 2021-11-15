---
title: "j2-1 Stoichiometry of lysine hydroxylations"
author: "Yoichiro Sugimoto"
date: "15 November, 2021"
vignette: >
  %\VignetteIndexEntry{Bioconductor style for PDF documents}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
   html_document:
     highlight: haddock
     toc: yes
     toc_depth: 2
     keep_md: yes
---


# Environment setup



```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)
library("Biostrings")

processors <- 8

j2.1.input.dir <- file.path("../../data/j2-PTM-stoichiometry/hydroxyK-stoichiometry")

results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")
j2.res.dir <- file.path(results.dir, "j2-PTM-stoichiometry")
j2.1.res.dir <- file.path(j2.res.dir, "j2-1-PTM-stoichiometry-K-only")

create.dirs(c(
    j2.res.dir,
    j2.1.res.dir
))
```


# Data import



```r
cells <- c(
    "HeLa_WT", "HeLa_J6KO",
    "HEK293", "MCF7",
    "HeLa_WT_J6PD", "HeLa_J6KO_J6PD", "HeLa_WT_J6PD_2"
)

all.pp.dt <- lapply(
    cells,
    function(x){
        dt <- fread(
            file.path(
                j2.1.input.dir,
                paste0(x, "__protein-peptides.csv")
            )
        )
        setnames(
            dt,
            old = c(
                "Protein Accession",
                grep("Area", colnames(dt), value = TRUE),
                "#Feature"
                ),
            new = c("Accession", "Area", "n_feature")
        )
        dt[, `:=`(
            cell = x,
            Peptide = gsub("(^[[:alpha:]]\\.|\\.[[:alpha:]]$)", "", Peptide)
        )]
        dt[, .(cell, Accession, Start, End, Unique, Peptide, Area, PTM, n_feature)]
    }
) %>% rbindlist

all.pp.dt[, cell := factor(cell, levels = cells)]

all.protein.feature.per.pos.dt <- file.path(
    j0.res.dir, "all_protein_feature_per_position.csv"
) %>% fread
```


# Data filtration

## Preliminary checks


```r
temp <- countPTM(all.pp.dt)
```

```
## [1] "All K related modifications"
## 
##                 K         K(+15.99) K(+15.99)(+42.01)         K(+42.01) 
##            105967              5687                 5               276 
## K(+42.01)(+15.99) K(+42.01)(+56.03) K(+42.01)(+72.02)         K(+56.03) 
##               122               189               130            244643 
## K(+56.03)(+42.01)         K(+72.02) K(+72.02)(+42.01) 
##                31             19490                 9 
## [1] "All R related modifications"
## 
##         R R(+42.01) 
##    216335        68 
## [1] "All PTMs"
## 
##         Acetylation (N-term) Acetylation (Protein N-term) 
##                         1923                          932 
##         Carbamidomethylation             Deamidation (NQ) 
##                         3303                        12292 
##                Oxidation (K)               Oxidation (MW) 
##                          241                        26153 
##      Oxidised Propionylation               Propionylation 
##                         2022                       109312
```

## Filtration and data processing

The following filtration and data processing will be performed:

1. Filter out non-unique peptide
2. Filter out peptides with Area == 0
3. Ignore acetylation
   - (e.g.) Treat K(+42.01) as K or K(+42.01)(+15.99) as K(+15.99)
4. Filter out peptides (1) with C-term propionylated K or dimethylated R and (2) without C-term K, hydroxyK, R or mono-methylated R
   - di-methylated R is not cleaved by trypsin 
5. Allow only 1 misclevage by trypsin 
   - At most one occurrence of K, K(+15.99), R, R(+14.02) in the peptide
   - Note that the following are allowed in the peptides
	 - C terminal K, R, KP, RP
	 - Propionylated K: K(+56.03), K(+72.02)
	 - Dimethyl R: R(+28.03)
6. The PTMs are assigned in the PTM column


```r
all.indexed.pp.dt <- filterPeptideData(all.pp.dt)
```

```
## [1] "Input peptides per cells"
## cell
##        HeLa_WT      HeLa_J6KO         HEK293           MCF7   HeLa_WT_J6PD 
##          24108          26136           5564           6982          73851 
## HeLa_J6KO_J6PD HeLa_WT_J6PD_2 
##          94827          24915 
## [1] "Peptides per cells after filteration 1"
## cell
##        HeLa_WT      HeLa_J6KO         HEK293           MCF7   HeLa_WT_J6PD 
##          22754          24032           4765           5969          66123 
## HeLa_J6KO_J6PD HeLa_WT_J6PD_2 
##          84717          22947 
## [1] "Peptides per cells after filtration 2"
## cell
##        HeLa_WT      HeLa_J6KO         HEK293           MCF7   HeLa_WT_J6PD 
##          18332          19075           4434           5414          59245 
## HeLa_J6KO_J6PD HeLa_WT_J6PD_2 
##          76529          20601 
## [1] "After the filtration of acetylation (filteration 3)"
## [1] "All K related modifications"
## 
##         K K(+15.99) K(+56.03) K(+72.02) 
##     81559      3241    193187     12487 
## [1] "All R related modifications"
## 
##      R 
## 170818 
## [1] "All PTMs"
## 
##         Acetylation (N-term) Acetylation (Protein N-term) 
##                         1474                          743 
##         Carbamidomethylation             Deamidation (NQ) 
##                         2461                         8454 
##                Oxidation (K)               Oxidation (MW) 
##                          172                        21514 
##      Oxidised Propionylation               Propionylation 
##                         1446                        89620 
## [1] "Peptides per cells after filtration 4"
## cell
##        HeLa_WT      HeLa_J6KO         HEK293           MCF7   HeLa_WT_J6PD 
##          17842          18657           4265           5267          57717 
## HeLa_J6KO_J6PD HeLa_WT_J6PD_2 
##          74667          19789 
## [1] "Peptides per cells after filtration 5"
## cell
##        HeLa_WT      HeLa_J6KO         HEK293           MCF7   HeLa_WT_J6PD 
##          16930          17507           4056           4946          56869 
## HeLa_J6KO_J6PD HeLa_WT_J6PD_2 
##          73559          19221 
## [1] "Filtration 5 (PTM assignment)"
## [1] "All K related modifications"
## 
##         K K(+15.99) K(+56.03) K(+72.02) 
##     82810       141    181006      1618 
## [1] "All R related modifications"
## 
##      R 
## 160854 
## [1] "All PTMs"
## 
##         Acetylation (N-term) Acetylation (Protein N-term) 
##                         1407                          710 
##         Carbamidomethylation             Deamidation (NQ) 
##                         2409                         8073 
##                Oxidation (K)               Oxidation (MW) 
##                          130                        20839 
##      Oxidised Propionylation               Propionylation 
##                         1303                        85896
```

# Individual amino acid position



```r
all.pos.dt <- lapply(
    c("K", "oxK", "R", "mR", "dmR", "all_aa"),
    indexToPosition,
    dt = all.indexed.pp.dt
) %>%
    rbindlist

print("The number of amino acids and PTM analysed / identified.")
```

```
## [1] "The number of amino acids and PTM analysed / identified."
```

```r
all.pos.dt[total_area > 0, table(cell, aa_type)]
```

```
##                 aa_type
## cell             all_aa      K    oxK      R
##   HeLa_WT        121453   9710     99   7722
##   HeLa_J6KO      127930   9977     48   8144
##   HEK293          24071   1891     56   1453
##   MCF7            30392   2492     47   1780
##   HeLa_WT_J6PD   412058  30130    327  30514
##   HeLa_J6KO_J6PD 516609  39168    272  36393
##   HeLa_WT_J6PD_2 143835  11085    224  11266
```



# Analysis of stoichiometry and other features

## K stoichiometry


```r
stoic.dt <- dcast(
    all.pos.dt[aa_type %in% c("K", "oxK")],
    formula = cell + Accession + position ~ aa_type,
    value.var = c("total_area", "total_n_feature"),
    fill = 0
)

stoic.dt[, `:=`(
    oxK_ratio = total_area_oxK / total_area_K
)]

all.protein.bs <- readAAStringSet("../../data/all_protein.fasta")
names(all.protein.bs) <- names(all.protein.bs) %>%
    {gsub("sp\\|", "", .)} %>%
    {str_split_fixed(., " ", n = 2)[, 1]}
uniprot.accession.flag <- names(all.protein.bs) %>%
    {str_split_fixed(., "\\|", n = 2)[, 2]} %>%
    {nchar(.) > 0}
all.protein.bs <- all.protein.bs[uniprot.accession.flag]

stoic.dt[
    , seq5 := BSgenome::getSeq(
                          x = all.protein.bs,
                          name = Accession
                        ) %>% as.character %>%
              substr(start = position - 5, stop = position + 5)
]

merge(
    all.protein.feature.per.pos.dt,
    stoic.dt,
    by = c("Accession", "position")
) %>%
    fwrite(file.path(j2.1.res.dir, "long_K_stoichiometry_data.csv"))

wide.stoic.dt <- dcast(
    stoic.dt,
    formula = Accession + position + seq5 ~ cell,
    value.var = c(
        "oxK_ratio",
        grep("^total_area", colnames(stoic.dt), value = TRUE),
        grep("^total_n_feature", colnames(stoic.dt), value = TRUE)
    )
)

merge(
    all.protein.feature.per.pos.dt,
    wide.stoic.dt,
    by = c("Accession", "position")
) %>%
    fwrite(
        file.path(j2.1.res.dir, "wide_K_stoichiometry_data.csv")
    )
```

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
##  date     2021-11-15                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package              * version  date       lib source        
##  Biobase                2.48.0   2020-04-27 [1] Bioconductor  
##  BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor  
##  BiocParallel           1.22.0   2020-04-27 [1] Bioconductor  
##  Biostrings           * 2.56.0   2020-04-27 [1] Bioconductor  
##  bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
##  BSgenome               1.56.0   2020-04-27 [1] Bioconductor  
##  cli                    3.0.1    2021-07-17 [1] CRAN (R 4.0.0)
##  colorspace             1.4-1    2019-03-18 [1] CRAN (R 4.0.0)
##  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
##  data.table           * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)
##  DelayedArray           0.14.0   2020-04-27 [1] Bioconductor  
##  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
##  dplyr                * 1.0.0    2020-05-29 [1] CRAN (R 4.0.0)
##  ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
##  evaluate               0.14     2019-05-28 [1] CRAN (R 4.0.0)
##  generics               0.0.2    2018-11-29 [1] CRAN (R 4.0.0)
##  GenomeInfoDb           1.24.0   2020-04-27 [1] Bioconductor  
##  GenomeInfoDbData       1.2.3    2021-06-10 [1] Bioconductor  
##  GenomicAlignments      1.24.0   2020-04-27 [1] Bioconductor  
##  GenomicRanges          1.40.0   2020-04-27 [1] Bioconductor  
##  ggplot2              * 3.3.1    2020-05-28 [1] CRAN (R 4.0.0)
##  glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
##  gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.0.0)
##  htmltools              0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
##  IRanges              * 2.22.1   2020-04-28 [1] Bioconductor  
##  khroma               * 1.3.0    2019-10-26 [1] CRAN (R 4.0.0)
##  knitr                * 1.28     2020-02-06 [1] CRAN (R 4.0.0)
##  lattice                0.20-41  2020-04-02 [1] CRAN (R 4.0.0)
##  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
##  magrittr             * 1.5      2014-11-22 [1] CRAN (R 4.0.0)
##  Matrix                 1.2-18   2019-11-27 [1] CRAN (R 4.0.0)
##  matrixStats            0.56.0   2020-03-13 [1] CRAN (R 4.0.0)
##  munsell                0.5.0    2018-06-12 [1] CRAN (R 4.0.0)
##  pillar                 1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
##  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
##  R6                     2.4.1    2019-11-12 [1] CRAN (R 4.0.0)
##  Rcpp                   1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
##  RCurl                  1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
##  rlang                  0.4.11   2021-04-30 [1] CRAN (R 4.0.0)
##  rmarkdown            * 2.2      2020-05-31 [1] CRAN (R 4.0.0)
##  Rsamtools              2.4.0    2020-04-27 [1] Bioconductor  
##  rtracklayer            1.48.0   2020-04-27 [1] Bioconductor  
##  S4Vectors            * 0.26.0   2020-04-27 [1] Bioconductor  
##  scales                 1.1.1    2020-05-11 [1] CRAN (R 4.0.0)
##  sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 4.0.5)
##  stringi                1.4.6    2020-02-17 [1] CRAN (R 4.0.0)
##  stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
##  SummarizedExperiment   1.18.1   2020-04-30 [1] Bioconductor  
##  tibble                 3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
##  tidyselect             1.1.0    2020-05-11 [1] CRAN (R 4.0.0)
##  vctrs                  0.3.1    2020-06-05 [1] CRAN (R 4.0.0)
##  withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
##  xfun                   0.14     2020-05-20 [1] CRAN (R 4.0.0)
##  XML                    3.99-0.3 2020-01-20 [1] CRAN (R 4.0.0)
##  XVector              * 0.28.0   2020-04-27 [1] Bioconductor  
##  yaml                   2.2.1    2020-02-01 [1] CRAN (R 4.0.0)
##  zlibbioc               1.34.0   2020-04-27 [1] Bioconductor  
## 
## [1] /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/hydroxylation_by_JMJD6/lib/R/library
```