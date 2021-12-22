j0-3. Summary of Protein feature analyses
================
Yoichiro Sugimoto
22 December, 2021

  - [Environment setup](#environment-setup)
  - [Merge all data](#merge-all-data)
  - [Export data](#export-data)
  - [Session information](#session-information)

# Environment setup

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

processors <- 8
```

``` r
results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")

create.dirs(c(
))
```

# Merge all data

``` r
all.iupred2.dt <- 
    file.path(j0.res.dir, "IUPRED2_per_protein_and_position.csv") %>%
    fread

k.ratio.dt <- 
    file.path(j0.res.dir, "K_ratio_per_protein_and_position.csv") %>%
    fread

all.biophysic.dt <- fread(file.path(j0.res.dir, "biophysical_property.csv"))
setnames(all.biophysic.dt, "Position", "position")

all.protein.feature.per.pos.dt <- merge(
    all.iupred2.dt,
    k.ratio.dt,
    by = c("Accession", "uniprot_id", "position"),
    all = TRUE
) %>%
    merge(
        all.biophysic.dt,
        by = c("Accession", "position"),
        all = TRUE
    )

brd.feature.per.pos.dt <- all.protein.feature.per.pos.dt[grepl("BRD[[:digit:]]", Accession)]
```

# Export data

``` r
fwrite(
    all.protein.feature.per.pos.dt,
    file.path(j0.res.dir, "all_protein_feature_per_position.csv")
)

fwrite(
    brd.feature.per.pos.dt,
    file.path(j0.res.dir, "BRDs_feature_per_position.csv")
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
    ##  date     2021-12-22                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date       lib source        
    ##  cli           3.0.1   2021-07-17 [1] CRAN (R 4.0.0)
    ##  colorspace    1.4-1   2019-03-18 [1] CRAN (R 4.0.0)
    ##  crayon        1.3.4   2017-09-16 [1] CRAN (R 4.0.0)
    ##  data.table  * 1.12.8  2019-12-09 [1] CRAN (R 4.0.0)
    ##  digest        0.6.25  2020-02-23 [1] CRAN (R 4.0.0)
    ##  dplyr       * 1.0.0   2020-05-29 [1] CRAN (R 4.0.0)
    ##  ellipsis      0.3.1   2020-05-15 [1] CRAN (R 4.0.0)
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 4.0.0)
    ##  generics      0.0.2   2018-11-29 [1] CRAN (R 4.0.0)
    ##  ggplot2     * 3.3.1   2020-05-28 [1] CRAN (R 4.0.0)
    ##  glue          1.4.1   2020-05-13 [1] CRAN (R 4.0.0)
    ##  gtable        0.3.0   2019-03-25 [1] CRAN (R 4.0.0)
    ##  htmltools     0.4.0   2019-10-04 [1] CRAN (R 4.0.0)
    ##  khroma      * 1.3.0   2019-10-26 [1] CRAN (R 4.0.0)
    ##  knitr       * 1.28    2020-02-06 [1] CRAN (R 4.0.0)
    ##  lifecycle     0.2.0   2020-03-06 [1] CRAN (R 4.0.0)
    ##  magrittr    * 1.5     2014-11-22 [1] CRAN (R 4.0.0)
    ##  munsell       0.5.0   2018-06-12 [1] CRAN (R 4.0.0)
    ##  pillar        1.4.4   2020-05-05 [1] CRAN (R 4.0.0)
    ##  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.0.0)
    ##  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.0.0)
    ##  R6            2.4.1   2019-11-12 [1] CRAN (R 4.0.0)
    ##  Rcpp          1.0.4.6 2020-04-09 [1] CRAN (R 4.0.0)
    ##  rlang         0.4.11  2021-04-30 [1] CRAN (R 4.0.0)
    ##  rmarkdown   * 2.2     2020-05-31 [1] CRAN (R 4.0.0)
    ##  scales        1.1.1   2020-05-11 [1] CRAN (R 4.0.0)
    ##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.0.5)
    ##  stringi       1.4.6   2020-02-17 [1] CRAN (R 4.0.0)
    ##  stringr     * 1.4.0   2019-02-10 [1] CRAN (R 4.0.0)
    ##  tibble        3.0.1   2020-04-20 [1] CRAN (R 4.0.0)
    ##  tidyselect    1.1.0   2020-05-11 [1] CRAN (R 4.0.0)
    ##  vctrs         0.3.1   2020-06-05 [1] CRAN (R 4.0.0)
    ##  withr         2.2.0   2020-04-20 [1] CRAN (R 4.0.0)
    ##  xfun          0.14    2020-05-20 [1] CRAN (R 4.0.0)
    ##  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.0.0)
    ## 
    ## [1] /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/hydroxylation_by_JMJD6/lib/R/library
