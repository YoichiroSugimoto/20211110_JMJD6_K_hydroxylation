j0-2 Protein feature extraction 2
================
Yoichiro Sugimoto
25 April, 2022

  - [Introduction](#introduction)
  - [Environment setup](#environment-setup)
  - [Data import](#data-import)
  - [Local hydropathy](#local-hydropathy)
  - [Local charge](#local-charge)
  - [Export data](#export-data)
  - [Session information](#session-information)

# Introduction

Additional biophysical properties of proteins will be analysed.

# Environment setup

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

library("Biostrings")
library("GenomicRanges")
library("idpr")

processors <- 8
```

``` r
results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")
ind.fa.dir <- file.path(j0.res.dir, "individual_fastas")

create.dirs(c(
))
```

# Data import

``` r
all.protein.bs <- readAAStringSet("../../data/all_protein.fasta")

names(all.protein.bs) <- names(all.protein.bs) %>%
    {gsub("sp\\|", "", .)} %>%
    {str_split_fixed(., " ", n = 2)[, 1]}

uniprot.accession.flag <- names(all.protein.bs) %>%
    {str_split_fixed(., "\\|", n = 2)[, 2]} %>%
    {nchar(.) > 0}

all.protein.bs <- all.protein.bs[uniprot.accession.flag]

## These 2 filtrations were peformed for this script only
all.protein.bs <- all.protein.bs[width(all.protein.bs) >= 9]

non.standard.aa.count <-lapply(
    all.protein.bs,
    function(x){
        str_count(as.character(x), "(U|O)")
    }
) %>% unlist

all.protein.bs <- all.protein.bs[non.standard.aa.count == 0]
```

# Local hydropathy

``` r
hydropathy.dt <- mcmapply(
    function(x, y){
        scaledHydropathyLocal(
            x,
            window = 11,
            plotResults = FALSE
        ) %>% data.table %>%
            {.[, Accession := y]}
    },
    all.protein.bs,
    names(all.protein.bs),
    SIMPLIFY = FALSE,
    mc.cores = processors
) %>%
    rbindlist
```

# Local charge

``` r
charge.dt <- mcmapply(
    function(x, y){
        chargeCalculationLocal(
            x,
            window = 11,
            plotResults = FALSE
        ) %>% data.table %>%
            {.[, Accession := y]}
    },
    all.protein.bs,
    names(all.protein.bs),
    SIMPLIFY = FALSE,
    mc.cores = processors
) %>%
    rbindlist
```

# Export data

``` r
all.biophysic.dt <- merge(
    hydropathy.dt[, .(Accession, Position, WindowHydropathy)],
    charge.dt[, .(Accession, Position, windowCharge, CenterResidue, Window)],
    by = c("Accession", "Position")
)

fwrite(
    all.biophysic.dt,
    file.path(
        j0.res.dir, "biophysical_property.csv"
    )
)

fwrite(
    charge.dt,
    file.path(
        j0.res.dir, "charge_per_position.csv"
    )
)
```

# Session information

``` r
sessioninfo::session_info()
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
    ##  date     2022-04-25                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package          * version  date       lib source        
    ##  BiocGenerics     * 0.34.0   2020-04-27 [1] Bioconductor  
    ##  Biostrings       * 2.56.0   2020-04-27 [1] Bioconductor  
    ##  bitops             1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
    ##  cli                3.0.1    2021-07-17 [1] CRAN (R 4.0.0)
    ##  colorspace         1.4-1    2019-03-18 [1] CRAN (R 4.0.0)
    ##  crayon             1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
    ##  data.table       * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)
    ##  digest             0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
    ##  dplyr            * 1.0.0    2020-05-29 [1] CRAN (R 4.0.0)
    ##  ellipsis           0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
    ##  evaluate           0.14     2019-05-28 [1] CRAN (R 4.0.0)
    ##  generics           0.0.2    2018-11-29 [1] CRAN (R 4.0.0)
    ##  GenomeInfoDb     * 1.24.0   2020-04-27 [1] Bioconductor  
    ##  GenomeInfoDbData   1.2.3    2021-06-10 [1] Bioconductor  
    ##  GenomicRanges    * 1.40.0   2020-04-27 [1] Bioconductor  
    ##  ggplot2          * 3.3.1    2020-05-28 [1] CRAN (R 4.0.0)
    ##  glue               1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
    ##  gtable             0.3.0    2019-03-25 [1] CRAN (R 4.0.0)
    ##  htmltools          0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
    ##  idpr             * 1.2.0    2021-05-19 [1] Bioconductor  
    ##  IRanges          * 2.22.1   2020-04-28 [1] Bioconductor  
    ##  khroma           * 1.3.0    2019-10-26 [1] CRAN (R 4.0.0)
    ##  knitr            * 1.28     2020-02-06 [1] CRAN (R 4.0.0)
    ##  lifecycle          0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
    ##  magrittr         * 1.5      2014-11-22 [1] CRAN (R 4.0.0)
    ##  munsell            0.5.0    2018-06-12 [1] CRAN (R 4.0.0)
    ##  pillar             1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
    ##  pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
    ##  purrr              0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
    ##  R6                 2.4.1    2019-11-12 [1] CRAN (R 4.0.0)
    ##  Rcpp               1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
    ##  RCurl              1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
    ##  rlang              0.4.11   2021-04-30 [1] CRAN (R 4.0.0)
    ##  rmarkdown        * 2.2      2020-05-31 [1] CRAN (R 4.0.0)
    ##  S4Vectors        * 0.26.0   2020-04-27 [1] Bioconductor  
    ##  scales             1.1.1    2020-05-11 [1] CRAN (R 4.0.0)
    ##  sessioninfo        1.1.1    2018-11-05 [1] CRAN (R 4.0.5)
    ##  stringi            1.4.6    2020-02-17 [1] CRAN (R 4.0.0)
    ##  stringr          * 1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
    ##  tibble             3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
    ##  tidyselect         1.1.0    2020-05-11 [1] CRAN (R 4.0.0)
    ##  vctrs              0.3.1    2020-06-05 [1] CRAN (R 4.0.0)
    ##  withr              2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
    ##  xfun               0.14     2020-05-20 [1] CRAN (R 4.0.0)
    ##  XVector          * 0.28.0   2020-04-27 [1] Bioconductor  
    ##  yaml               2.2.1    2020-02-01 [1] CRAN (R 4.0.0)
    ##  zlibbioc           1.34.0   2020-04-27 [1] Bioconductor  
    ## 
    ## [1] /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/hydroxylation_by_JMJD6/lib/R/library
