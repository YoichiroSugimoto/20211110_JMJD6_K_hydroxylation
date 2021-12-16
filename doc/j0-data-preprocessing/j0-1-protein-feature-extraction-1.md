j0-1. Protein feature extraction (1/2)
================
Yoichiro Sugimoto
16 December, 2021

  - [Overview](#overview)
  - [Environment setup](#environment-setup)
  - [Data import](#data-import)
  - [IUPRED2](#iupred2)
  - [K ratio](#k-ratio)
  - [Session information](#session-information)

# Overview

In this script, disorderedness and enrichment of lysine residues in
protein sequences will be examined.

# Environment setup

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

library("Biostrings")
library("GenomicRanges")

processors <- 8
```

``` r
results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")

ind.fa.dir <- file.path(j0.res.dir, "individual_fastas")
iupred2.dir <- normalizePath(file.path(j0.res.dir, "IUPred2_score"))
```

    ## Warning in normalizePath(file.path(j0.res.dir, "IUPred2_score")):
    ## path[1]="/camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/projects/
    ## 20211110_JMJD6_K_hydroxylation/results/j0-data-preprocessing/IUPred2_score": No
    ## such file or directory

``` r
create.dirs(c(
    results.dir,
    j0.res.dir,
    ind.fa.dir,
    iupred2.dir
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
```

# IUPRED2

``` r
temp <- mclapply(
    1:length(all.protein.bs),
    function(x){
        writeXStringSet(
            all.protein.bs[x],
            file.path(
                ind.fa.dir,
                paste0(
                    str_split_fixed(
                        names(all.protein.bs)[x], "\\|", 2
                    )[, 1], ".fa")
            )
        )
    },
    mc.cores = processors
)

temp <- mclapply(
    1:length(all.protein.bs),
    function(x){
        paste(
            "python",
            file.path("/camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/source_files/iupred2a/iupred2a.py"),
            file.path(
                ind.fa.dir,
                paste0(
                    str_split_fixed(
                        names(all.protein.bs)[x], "\\|", 2
                    )[, 1], ".fa")
            ),
            "long", ">",
            file.path(
                iupred2.dir,
                paste0(
                    str_split_fixed(
                        names(all.protein.bs)[x], "\\|", 2
                    )[, 1], ".txt"
                )
            )
        ) %>%
            system
    },
    mc.cores = processors
)

all.iupred2.dt <- mcmapply(
    function(x, y){
        dt <- fread(file.path(iupred2.dir, paste0(x, ".txt")))
        setnames(
            dt,
            old = c("# POS", "RES", "IUPRED2"),
            new = c("position", "residue", "IUPRED2")
        )
        dt[, `:=`(
            Accession = y,
            uniprot_id = x
        )]
    },
    str_split_fixed(names(all.protein.bs), "\\|", n = 2)[, 1],
    names(all.protein.bs),
    SIMPLIFY = FALSE,
    mc.cores = processors
) %>%
    rbindlist

fwrite(
    all.iupred2.dt,
    file.path(j0.res.dir, "IUPRED2_per_protein_and_position.csv")
)


max.iupred2.dt <-
    copy(all.iupred2.dt)[, list(
            IUPRED2_max = max(IUPRED2),
            IUPRED2_mean = mean(IUPRED2)
        ), by = list(Accession, uniprot_id)]

fwrite(
    max.iupred2.dt,
    file.path(j0.res.dir, "Max_IUPRED2_per_protein.csv")
)

temp <- do.call(file.remove, list(list.files(ind.fa.dir, full.names = TRUE)))
temp <- do.call(file.remove, list(list.files(iupred2.dir, full.names = TRUE)))
```

# K ratio

``` r
calculateLettersIn10merRatio <- function(sl.aa.letters, all.protein.bs, window.size = 10){
    
    letter.pos.list <- mclapply(
        all.protein.bs,
        letterFrequencyInSlidingView,
        view.width = 1, letters = sl.aa.letters, as.prob = TRUE,
        mc.cores = processors
    )

    ratio.dt <- mcmapply(
        function(x, y){
            data.table(
                position = 1:nrow(x),
                letter_position = rowSums(x[, sl.aa.letters, drop = FALSE]),
                letter_ratio = zoo::rollapply(
                                        rowSums(x[, sl.aa.letters, drop = FALSE]),
                                        window.size,
                                        FUN = function(y){sum(y) / window.size},
                                        partial = TRUE, align = "left"
                                    ),
                Accession = y
            )
        },
        letter.pos.list,
        names(letter.pos.list),
        SIMPLIFY = FALSE,
        mc.cores = processors
    ) %>%
        rbindlist

    ratio.dt[, uniprot_id := str_split_fixed(Accession, "\\|", n = 2)[, 1]]

    calcMaxLetter <- function(accession, dt){
        dt <- copy(dt)[Accession == accession]
        dt[, `:=`(
            letter_ratio_score =
                data.table::shift(letter_ratio, window.size - 1, "lag", fill = 0) %>%
                {zoo::rollapply(., window.size, FUN = max, partial = TRUE, align = "left")}
        )]
        return(dt)
    }

    letter.ratio.dt <- mclapply(
        ratio.dt[, unique(Accession)],
        calcMaxLetter,
        dt = ratio.dt,
        mc.cores = processors
    ) %>%
        rbindlist

    max.letter.ratio.dt <-
        copy(letter.ratio.dt)[, list(
                letter_ratio_max = max(letter_ratio)
            ), by = list(Accession, uniprot_id)]

    setnames(
        letter.ratio.dt,
        old = colnames(letter.ratio.dt),
        new = gsub("letter", paste(sl.aa.letters, collapse = ""), colnames(letter.ratio.dt))
    )

    setnames(
        max.letter.ratio.dt,
        old = colnames(max.letter.ratio.dt),
        new = gsub("letter", paste(sl.aa.letters, collapse = ""), colnames(max.letter.ratio.dt))
    )

    fwrite(
        letter.ratio.dt,
        file.path(j0.res.dir, paste0(
                                   paste(sl.aa.letters, collapse = ""),
                                   "_ratio_per_protein_and_position.csv")
                  )
    )

    fwrite(
        max.letter.ratio.dt,
        file.path(
            j0.res.dir,
            paste0(
                "Max_",
                paste(sl.aa.letters, collapse = ""),
                "_ratio_per_protein.csv"
            )
        )
    )

    return()
}
```

``` r
temp <- calculateLettersIn10merRatio(
    sl.aa.letters = "K",
    all.protein.bs = all.protein.bs,
    window.size = 10
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
    ##  date     2021-12-16                  
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
