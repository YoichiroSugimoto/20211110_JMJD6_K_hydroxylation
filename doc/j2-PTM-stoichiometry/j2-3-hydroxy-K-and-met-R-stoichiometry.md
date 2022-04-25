j2-3 Analysis of the stoichiometry of lysine hydroxylation and arginine
methylation
================
Yoichiro Sugimoto
22 April, 2022

  - [Environment setup](#environment-setup)
  - [Data import](#data-import)
  - [Data filtration](#data-filtration)
      - [Preliminary checks](#preliminary-checks)
      - [Filtration and data
        processing](#filtration-and-data-processing)
  - [Analysis of stoichiometry and other
    features](#analysis-of-stoichiometry-and-other-features)
      - [K stoichiometry](#k-stoichiometry)
      - [R stoichiometry](#r-stoichiometry)
  - [Session information](#session-information)

# Environment setup

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)
temp <- sapply(list.files("./functions", full.names = TRUE), source)
library("Biostrings")

processors <- 8

j2.input.dir <- file.path("../../data/j2-PTM-stoichiometry")
j2.3.input.dir <- file.path(j2.input.dir, "hydroxyK-methylR-stoichiometry")

results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")
j2.res.dir <- file.path(results.dir, "j2-PTM-stoichiometry")
j2.3.res.dir <- file.path(j2.res.dir, "j2-3-PTM-stoichiometry-K-and-R")

create.dirs(c(
    j2.res.dir,
    j2.3.res.dir
))
```

# Data import

``` r
data.source <- c(
    "HeLa_WT_JQ1",
    "HeLa_JMJD6KO_JQ1",
    "HeLa_WT_J6pep",
    "HeLa_JMJD6KO_J6pep"
)

all.pp.dt<- lapply(
    1:length(data.source),
    function(x){
        dt <- fread(
            paste0(
                j2.3.input.dir,
                "/",
                data.source[x], "_Trypsin_oxKandMetR__protein-peptides.csv"
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
            data_source = data.source[x],
            Peptide = gsub("(^[[:alpha:]]\\.|\\.[[:alpha:]]$)", "", Peptide)
        )]
        dt[, .(data_source, Accession, Start, End, Unique, Peptide, Area, PTM, n_feature)]
    }
) %>% rbindlist

all.pp.dt[, data_source := factor(data_source, levels = data.source)]

print("The number of peptides per data source")
```

    ## [1] "The number of peptides per data source"

``` r
all.pp.dt[, .N, by = data_source]
```

    ##           data_source     N
    ## 1:        HeLa_WT_JQ1 24502
    ## 2:   HeLa_JMJD6KO_JQ1 26774
    ## 3:      HeLa_WT_J6pep 75211
    ## 4: HeLa_JMJD6KO_J6pep 96360

``` r
all.protein.feature.per.pos.dt <- file.path(
    j0.res.dir, "all_protein_feature_per_position.csv"
) %>% fread
```

# Data filtration

## Preliminary checks

``` r
temp <- countPTM(all.pp.dt)
```

    ## [1] "All K related modifications"
    ## 
    ##                 K         K(+15.99) K(+15.99)(+42.01)         K(+42.01) 
    ##             90621              4000                 1               306 
    ## K(+42.01)(+15.99) K(+42.01)(+56.03) K(+42.01)(+72.02)         K(+56.03) 
    ##               123               197               131            221758 
    ## K(+56.03)(+42.01)         K(+72.02) K(+72.02)(+42.01) 
    ##                19             16666                 5 
    ## [1] "All R related modifications"
    ## 
    ##                 R         R(+14.02)         R(+28.03)         R(+42.01) 
    ##            188757              1713              1580                64 
    ## R(+42.01)(+14.02) R(+42.01)(+28.03) 
    ##                74                17 
    ## [1] "All PTMs"
    ## 
    ##         Acetylation (N-term) Acetylation (Protein N-term) 
    ##                         2011                          646 
    ##         Carbamidomethylation             Deamidation (NQ) 
    ##                         3339                        10575 
    ##             Dimethylation(R)               Methylation(R) 
    ##                          294                          238 
    ##                Oxidation (K)               Oxidation (MW) 
    ##                          139                        22938 
    ##      Oxidised Propionylation               Propionylation 
    ##                         1572                        99375

## Filtration and data processing

The following filtration and data processing will be performed:

1.  Filter out non-unique peptide
2.  Filter out peptides with Area == 0
3.  Ignore acetylation
      - (e.g.) Treat K(+42.01) as K or K(+42.01)(+15.99) as K(+15.99)
4.  Filter out peptides (1) with C-term propionylated K or dimethylated
    R and (2) without C-term K, hydroxyK, R or mono-methylated R
      - di-methylated R is not cleaved by trypsin
5.  Allow only 1 misclevage by trypsin
      - At most one occurrence of K, K(+15.99), R, R(+14.02) in the
        peptide
      - Note that the following are allowed in the peptides
          - C terminal K and R, KP and RP
          - Propionylated K: K(+56.03), K(+72.02)
          - Dimethyl R: R(+28.03)
6.  The PTMs are assigned in the PTM column
7.  (additional filter) K does not have M or W within 2 amino acids to
    be assigned as hydroxy K

<!-- end list -->

``` r
all.indexed.pp.dt <- filterPeptideData(all.pp.dt)
```

    ## [1] "Input peptides per data_sources"
    ## data_source
    ##        HeLa_WT_JQ1   HeLa_JMJD6KO_JQ1      HeLa_WT_J6pep HeLa_JMJD6KO_J6pep 
    ##              24502              26774              75211              96360 
    ## [1] "Peptides per data_sources after filteration 1"
    ## data_source
    ##        HeLa_WT_JQ1   HeLa_JMJD6KO_JQ1      HeLa_WT_J6pep HeLa_JMJD6KO_J6pep 
    ##              23020              24538              67124              85954 
    ## [1] "Peptides per data_sources after filtration 2"
    ## data_source
    ##        HeLa_WT_JQ1   HeLa_JMJD6KO_JQ1      HeLa_WT_J6pep HeLa_JMJD6KO_J6pep 
    ##              18516              19349              60023              77538 
    ## [1] "After the filtration of acetylation (filteration 3)"
    ## [1] "All K related modifications"
    ## 
    ##         K K(+15.99) K(+56.03) K(+72.02) 
    ##     69195      2164    172614     10317 
    ## [1] "All R related modifications"
    ## 
    ##         R R(+14.02) R(+28.03) 
    ##    147749      1038       975 
    ## [1] "All PTMs"
    ## 
    ##         Acetylation (N-term) Acetylation (Protein N-term) 
    ##                         1513                          482 
    ##         Carbamidomethylation             Deamidation (NQ) 
    ##                         2475                         7230 
    ##             Dimethylation(R)               Methylation(R) 
    ##                          176                          157 
    ##                Oxidation (K)               Oxidation (MW) 
    ##                           93                        18729 
    ##      Oxidised Propionylation               Propionylation 
    ##                         1034                        80462 
    ## [1] "Peptides per data_sources after filtration 4"
    ## data_source
    ##        HeLa_WT_JQ1   HeLa_JMJD6KO_JQ1      HeLa_WT_J6pep HeLa_JMJD6KO_J6pep 
    ##              18015              18908              58412              75597 
    ## [1] "Peptides per data_sources after filtration 5"
    ## data_source
    ##        HeLa_WT_JQ1   HeLa_JMJD6KO_JQ1      HeLa_WT_J6pep HeLa_JMJD6KO_J6pep 
    ##              17074              17695              57449              74373 
    ## [1] "Filtration 5 (PTM assignment)"
    ## [1] "All K related modifications"
    ## 
    ##         K K(+15.99) K(+56.03) K(+72.02) 
    ##     70443        76    161973      1098 
    ## [1] "All R related modifications"
    ## 
    ##         R R(+14.02) R(+28.03) 
    ##    140338       147       107 
    ## [1] "All PTMs"
    ## 
    ##         Acetylation (N-term) Acetylation (Protein N-term) 
    ##                         1439                          458 
    ##         Carbamidomethylation             Deamidation (NQ) 
    ##                         2422                         6927 
    ##             Dimethylation(R)               Methylation(R) 
    ##                          107                          145 
    ##                Oxidation (K)               Oxidation (MW) 
    ##                           70                        18130 
    ##      Oxidised Propionylation               Propionylation 
    ##                          941                        77054

``` r
all.pos.dt <- lapply(
    c("K", "oxK", "R", "mR", "dmR", "all_aa"),
    indexToPosition,
    dt = all.indexed.pp.dt
) %>%
    rbindlist
```

# Analysis of stoichiometry and other features

## K stoichiometry

``` r
k.stoic.dt <- dcast(
    all.pos.dt[aa_type %in% c("K", "oxK")],
    formula = data_source + Accession + position ~ aa_type,
    value.var = c("total_area", "total_n_feature"),
    fill = 0
)

k.stoic.dt[, `:=`(
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

a.range <- 5

k.stoic.dt[
    , seq5 := BSgenome::getSeq(
                          x = all.protein.bs,
                          name = Accession
                        ) %>% as.character %>%
              substr(start = position - a.range, stop = position + a.range)
]

## print("Before filtration 7")
## k.stoic.dt[, table(Hydroxy_K = oxK_ratio > 0, data_source)]

## MW filter
k.stoic.dt[, `:=`(
    MW_within_1 = substr(
        seq5,
        start = min(a.range, position - 1),
        stop = min(a.range, position - 1) + 2
    ) %>%
        str_count(pattern = "(M|W)") > 0,
    MW_within_2 = substr(
        seq5,
        start = min(a.range, position - 2),
        stop = min(a.range, position - 2) + 4
    ) %>%
        str_count(pattern = "(M|W)") > 0    
), by = seq_len(nrow(k.stoic.dt))]

## Import curated data
curated.hydroxK.site.dt <- fread(
    file.path(j2.input.dir, "curated_JMJD6_substrate_sites.csv") 
)

curated.hydroxK.site.dt[, `:=`(
    curated_oxK_site = TRUE
)]

curated.hydroxK.site.dt <- curated.hydroxK.site.dt[, .(
    Accession, position, curated_oxK_site, screen
)]

all.protein.feature.per.pos.dt <- merge(
    curated.hydroxK.site.dt,
    all.protein.feature.per.pos.dt,
    all.y = TRUE,
    by = c("Accession", "position")
)

all.protein.feature.per.pos.dt[, `:=`(
    curated_oxK_site = case_when(
        is.na(curated_oxK_site) ~ FALSE,
        TRUE ~ curated_oxK_site
    )
)]

merge(
    all.protein.feature.per.pos.dt,
    k.stoic.dt,
    by = c("Accession", "position")
) %>%
    fwrite(file.path(j2.3.res.dir, "long_K_stoichiometry_data.csv"))

wide.k.stoic.dt <- dcast(
    k.stoic.dt,
    formula = Accession + position + seq5 + MW_within_1 + MW_within_2 ~ data_source,
    value.var = c(
        "oxK_ratio",
        grep("^total_area", colnames(k.stoic.dt), value = TRUE),
        grep("^total_n_feature", colnames(k.stoic.dt), value = TRUE)
    )
)

merge(
    all.protein.feature.per.pos.dt,
    wide.k.stoic.dt,
    by = c("Accession", "position")
) %>%
    fwrite(
        file.path(j2.3.res.dir, "wide_K_stoichiometry_data.csv")
    )
```

## R stoichiometry

``` r
stoic.dt <- dcast(
    all.pos.dt[aa_type %in% c("R", "mR", "dmR")],
    formula = data_source + Accession + position ~ aa_type,
    value.var = c("total_area", "total_n_feature"),
    fill = 0
)

stoic.dt[, `:=`(
    mR_ratio = total_area_mR / total_area_R,
    dmR_ratio = total_area_dmR / total_area_R
)]

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
    fwrite(file.path(j2.3.res.dir, "long_R_stoichiometry_data.csv"))

wide.stoic.dt <- dcast(
    stoic.dt,
    formula = Accession + position + seq5 ~ data_source,
    value.var = c(
        "mR_ratio", "dmR_ratio",
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
        file.path(j2.3.res.dir, "wide_R_stoichiometry_data.csv")
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
    ##  date     2022-04-22                  
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
