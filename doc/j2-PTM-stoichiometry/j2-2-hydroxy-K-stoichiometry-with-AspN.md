Analysis of lysine hydroxylation stoichiometry for Asp-N data
================
Yoichiro Sugimoto
22 December, 2021

  - [Environment setup](#environment-setup)
  - [Data import](#data-import)
  - [Individual amino acid position](#individual-amino-acid-position)
  - [Analysis of stoichiometry and other
    features](#analysis-of-stoichiometry-and-other-features)
  - [Session information](#session-information)

The following values were manually calculated to double check the
outputs.

BRD4 DOX+ manually calculated stoichiometry Position All hydroxy\_K %
519 2.91E+07 0 0% 535 3.78E+08 1.23E+08 33% 537 3.78E+08 2.87E+08 76%
538 3.78E+08 2.70E+08 71% 539 3.78E+08 1.21E+08 32% 541 3.78E+08
2.25E+08 59%

# Environment setup

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)
library("Biostrings")
processors <- 8

j2.2.input.dir <- file.path("../../data/j2-PTM-stoichiometry/hydroxyK-stoichiometry_AspN")

results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")
j2.res.dir <- file.path(results.dir, "j2-PTM-stoichiometry")
j2.2.res.dir <- file.path(j2.res.dir, "j2-2-PTM-stoichiometry-K-only-AspN")

create.dirs(c(
    j2.res.dir,
    j2.2.res.dir
))
```

# Data import

``` r
pp.dt <- fread(
    file.path(
        j2.2.input.dir,
        "HeLa_WT_BRD4IP_AspN_oxK__protein-peptides.csv"
    )
)

original.data.colnames <- c(
    "MC48; No Dox", "MC50; Plus Dox; 21% O2", "MC52; Plus Dox; 1% O2", "MC54; Plus Dox; 0.1% O2"
)

new.data.colnames <- c(
    "Dox_minus__Ox21", "Dox_plus_Ox21", "Dox_plus_Ox1", "Dox_plus_Ox01"
)

setnames(
    pp.dt,
    old = c(
        "Protein Accession",
        paste("Area", original.data.colnames),
        paste("#Feature", original.data.colnames)
    ),
    new = c(
        "Accession",
        paste0("Area__", new.data.colnames),
        paste0("n_feature__", new.data.colnames)
    )
)

pp.dt[, `:=`(
    Peptide = gsub("(^[[:alpha:]]\\.|\\.[[:alpha:]]$)", "", Peptide),
    PTMs = str_split(PTM, "; ")
)]

pp.dt <- pp.dt[, c(
    "Accession", "Start", "End", "Unique", "Peptide", 
    paste0("Area__", new.data.colnames), paste0("n_feature__", new.data.colnames),
    "PTMs"
), with = FALSE]

## Only analyse unique peptides and peptides from BRD4
pp.dt <- pp.dt[
    Unique == "Y" &
    Accession == "O60885|BRD4_HUMAN" ##&
    ## End == 541
]

area.pp.dt <- melt(
    pp.dt,
    id.vars = c("Accession", "Peptide", "Start", "End", "Unique", "PTMs"),
    measure.vars = paste0("Area__", new.data.colnames),
    value.name = "Area"
) %>%
    {.[, data_source := str_split_fixed(variable, "__", n = 2)[, 2]]}

n.feature.pp.dt <- melt(
    pp.dt,
    id.vars = c("Accession", "Peptide", "Start", "End"),
    measure.vars = paste0("n_feature__", new.data.colnames),
    value.name = "n_feature"
) %>%
    {.[, data_source := str_split_fixed(variable, "__", n = 2)[, 2]]}


all.pp.dt <- merge(
    area.pp.dt,
    n.feature.pp.dt,
    by = c("Accession", "Peptide", "Start", "End", "data_source")
)

all.pp.dt <- all.pp.dt[
    !is.na(Area) & Area > 0
]


all.pp.dt[
  , corrected_peptide := gsub("\\(\\+\\.98\\)", "", Peptide) %>% # Remove (+.98)
        {gsub("\\(\\+42\\.01\\)", "", .)} # Remove Acetylation from modification
]

all.pp.dt[
  , corrected_peptide := case_when(
        !any(PTMs[[1]] == "Oxidation (K)") ~
            gsub("K\\(\\+15\\.99\\)","K", corrected_peptide),
        TRUE ~ corrected_peptide
    ) ## %>% {case_when(
      ##          !any(PTMs[[1]] == "Oxidised Propionylation") ~
      ##              gsub("K\\(\\+72.02\\)","K", .),
      ##          TRUE ~ .
    ##      )}
   ,
    by = seq_len(nrow(all.pp.dt))
]

all.pp.dt[, `:=`(
    raw_peptide = gsub("\\s*\\([^\\)]+\\)","", Peptide),
    k_info_peptide =
        gsub(
            names(AMINO_ACID_CODE) %>%
            {.[. != "K"]} %>%
            {paste(., collapse = "|")} %>%
            {paste0("(", ., ")")}, 
            "A", corrected_peptide # Non K amino acids will be replaced with A
        ) %>%
        {gsub("K\\(\\+56\\.03\\)", "K", .)} %>% 
        {gsub("(K\\(\\+15\\.99\\)|K\\(\\+72\\.02\\))", "H", .)} %>% # H: hydroxylared K
        {gsub("\\s*\\([^\\)]+\\)","", .)}
)]

all.pp.dt[, `:=`(
    K_index = str_locate_all(k_info_peptide, "(K|H)") %>%
        {lapply(., function(x) x[,2])},
    oxK_index = str_locate_all(k_info_peptide, "H") %>%
        {lapply(., function(x) x[,2])},
    all_aa_index = str_locate_all(k_info_peptide, "(K|H|A)") %>%
        {lapply(., function(x) x[,2])}
)]
```

# Individual amino acid position

``` r
all.pos.dt <- lapply(
    c("K", "oxK", "all_aa"),
    indexToPosition,
    dt = all.pp.dt
) %>%
    rbindlist
```

# Analysis of stoichiometry and other features

``` r
stoic.dt <- dcast(
    all.pos.dt[aa_type %in% c("K", "oxK")],
    formula = data_source + Accession + position ~ aa_type,
    value.var = c("total_area", "total_n_feature"),
    fill = 0
)

stoic.dt[, `:=`(
    ox_ratio = total_area_oxK / total_area_K
)]


wide.dt <- dcast(
    stoic.dt,
    Accession + position ~ data_source,
    value.var = c("ox_ratio", "total_n_feature_K", "total_n_feature_oxK")
)


fwrite(
    wide.dt,
    file.path(j2.2.res.dir, "hydroxylation_stoichiometry_for_AspN_data.csv")
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
    ##  package      * version date       lib source        
    ##  BiocGenerics * 0.34.0  2020-04-27 [1] Bioconductor  
    ##  Biostrings   * 2.56.0  2020-04-27 [1] Bioconductor  
    ##  cli            3.0.1   2021-07-17 [1] CRAN (R 4.0.0)
    ##  colorspace     1.4-1   2019-03-18 [1] CRAN (R 4.0.0)
    ##  crayon         1.3.4   2017-09-16 [1] CRAN (R 4.0.0)
    ##  data.table   * 1.12.8  2019-12-09 [1] CRAN (R 4.0.0)
    ##  digest         0.6.25  2020-02-23 [1] CRAN (R 4.0.0)
    ##  dplyr        * 1.0.0   2020-05-29 [1] CRAN (R 4.0.0)
    ##  ellipsis       0.3.1   2020-05-15 [1] CRAN (R 4.0.0)
    ##  evaluate       0.14    2019-05-28 [1] CRAN (R 4.0.0)
    ##  generics       0.0.2   2018-11-29 [1] CRAN (R 4.0.0)
    ##  ggplot2      * 3.3.1   2020-05-28 [1] CRAN (R 4.0.0)
    ##  glue           1.4.1   2020-05-13 [1] CRAN (R 4.0.0)
    ##  gtable         0.3.0   2019-03-25 [1] CRAN (R 4.0.0)
    ##  htmltools      0.4.0   2019-10-04 [1] CRAN (R 4.0.0)
    ##  IRanges      * 2.22.1  2020-04-28 [1] Bioconductor  
    ##  khroma       * 1.3.0   2019-10-26 [1] CRAN (R 4.0.0)
    ##  knitr        * 1.28    2020-02-06 [1] CRAN (R 4.0.0)
    ##  lifecycle      0.2.0   2020-03-06 [1] CRAN (R 4.0.0)
    ##  magrittr     * 1.5     2014-11-22 [1] CRAN (R 4.0.0)
    ##  munsell        0.5.0   2018-06-12 [1] CRAN (R 4.0.0)
    ##  pillar         1.4.4   2020-05-05 [1] CRAN (R 4.0.0)
    ##  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.0.0)
    ##  purrr          0.3.4   2020-04-17 [1] CRAN (R 4.0.0)
    ##  R6             2.4.1   2019-11-12 [1] CRAN (R 4.0.0)
    ##  Rcpp           1.0.4.6 2020-04-09 [1] CRAN (R 4.0.0)
    ##  rlang          0.4.11  2021-04-30 [1] CRAN (R 4.0.0)
    ##  rmarkdown    * 2.2     2020-05-31 [1] CRAN (R 4.0.0)
    ##  S4Vectors    * 0.26.0  2020-04-27 [1] Bioconductor  
    ##  scales         1.1.1   2020-05-11 [1] CRAN (R 4.0.0)
    ##  sessioninfo    1.1.1   2018-11-05 [1] CRAN (R 4.0.5)
    ##  stringi        1.4.6   2020-02-17 [1] CRAN (R 4.0.0)
    ##  stringr      * 1.4.0   2019-02-10 [1] CRAN (R 4.0.0)
    ##  tibble         3.0.1   2020-04-20 [1] CRAN (R 4.0.0)
    ##  tidyselect     1.1.0   2020-05-11 [1] CRAN (R 4.0.0)
    ##  vctrs          0.3.1   2020-06-05 [1] CRAN (R 4.0.0)
    ##  withr          2.2.0   2020-04-20 [1] CRAN (R 4.0.0)
    ##  xfun           0.14    2020-05-20 [1] CRAN (R 4.0.0)
    ##  XVector      * 0.28.0  2020-04-27 [1] Bioconductor  
    ##  yaml           2.2.1   2020-02-01 [1] CRAN (R 4.0.0)
    ##  zlibbioc       1.34.0  2020-04-27 [1] Bioconductor  
    ## 
    ## [1] /camp/lab/ratcliffep/home/users/sugimoy/CAMP_HPC/software/miniconda3_20200606/envs/hydroxylation_by_JMJD6/lib/R/library
