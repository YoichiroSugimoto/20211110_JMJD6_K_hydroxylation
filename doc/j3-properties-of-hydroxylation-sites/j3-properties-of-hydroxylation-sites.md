Analysis of lysine hydroxylation stoichiometry and protein features
================
Yoichiro Sugimoto
25 April, 2022

  - [Package import](#package-import)
  - [Data import](#data-import)
      - [Annotation](#annotation)
      - [Stoichiometry data](#stoichiometry-data)
  - [Analysis of the properties of hydroxylation
    sites](#analysis-of-the-properties-of-hydroxylation-sites)
  - [Analysis of sequence feature around hydroxylation
    site](#analysis-of-sequence-feature-around-hydroxylation-site)
      - [Amino acid enrichment around the hydroxylation
        sites](#amino-acid-enrichment-around-the-hydroxylation-sites)
      - [Biophysical properties](#biophysical-properties)
  - [Session information](#session-information)

# Package import

``` r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

library("Biostrings")
library("GenomicRanges")
library("BSgenome")
library("readxl")

processors <- 8

set.seed(1)
```

``` r
j3.input.dir <- file.path("../../data/j3-properties-of-hydroxylation-sites")

results.dir <- normalizePath(file.path("../../results"))
j0.res.dir <- file.path(results.dir, "j0-data-preprocessing")
j2.res.dir <- file.path(results.dir, "j2-PTM-stoichiometry")
j2.1.res.dir <- file.path(j2.res.dir, "j2-1-PTM-stoichiometry-K-only")

create.dirs(c(
))
```

# Data import

## Annotation

``` r
all.protein.bs <- readAAStringSet("../../data/all_protein.fasta")
names(all.protein.bs) <- names(all.protein.bs) %>%
    {gsub("sp\\|", "", .)} %>%
    {str_split_fixed(., " ", n = 2)[, 1]}
uniprot.accession.flag <- names(all.protein.bs) %>%
    {str_split_fixed(., "\\|", n = 2)[, 2]} %>%
    {nchar(.) > 0}
all.protein.bs <- all.protein.bs[uniprot.accession.flag]

all.protein.feature.per.pos.dt <- file.path(
    j0.res.dir, "all_protein_feature_per_position.csv"
) %>% fread
```

## Stoichiometry data

``` r
stoichiometry.dt <- fread(
    file.path(
        j2.1.res.dir,
        "long_K_stoichiometry_data.csv"
    )
)

stoichiometry.dt[
  , curated_oxK_site := case_when(
        curated_oxK_site == TRUE ~ "JMJD6_substrate",
        TRUE ~ "Others"
    )
]

stoichiometry.dt[, `:=`(
    stoichiometry_available = case_when(
        is.na(oxK_ratio) ~ FALSE,
        TRUE ~ TRUE
    ),
    Accession_position = paste(Accession, position, sep = "_")
)]

## Sanity check
j6pep.substrate.dt <- stoichiometry.dt[
    stoichiometry_available == TRUE &
    curated_oxK_site == "JMJD6_substrate" &
    data_source == "HeLa_WT_J6pep" &
    total_n_feature_K > 1 &
    grepl("J6pep", screen),
    .(Accession, screen, position, total_n_feature_K, oxK_ratio)
][order(oxK_ratio, decreasing = TRUE)]

j6pep.substrate.dt[
  , accession_position := str_split_fixed(Accession, "\\|", n = 2)[, 2] %>%
        {gsub("_HUMAN", "", .)} %>%
        {paste0(., " (", position, ")")} %>%
        {factor(., levels = .)}
]


reportOxKSiteStats <- function(stoichiometry.dt){
    print("Total unique JMJD6 substrate proteins and sites:")
    stoichiometry.dt[, list(
        unique_protein_N = uniqueN(Accession),
        unique_site_N = uniqueN(Accession_position)
    ),
    by = list(
        curated_oxK_site
    )
    ][order(curated_oxK_site)] %>%
        print

    print("by the availability of stocihiometry data")
    stoichiometry.dt[, list(
        unique_protein_N = uniqueN(Accession),
        unique_site_N = uniqueN(Accession_position)
    ),
    by = list(
        curated_oxK_site,
        stoichiometry_available
    )
    ][order(curated_oxK_site, -stoichiometry_available)] %>%
        print

    print("by the data source")
    stoichiometry.dt[, list(
        unique_protein_N = uniqueN(Accession),
        unique_site_N = uniqueN(Accession_position)
    ),
    by = list(
        curated_oxK_site,
        stoichiometry_available,
        data_source
    )
    ][order(curated_oxK_site, data_source)] %>%
        print

    return()
}

temp <- reportOxKSiteStats(stoichiometry.dt)
```

    ## [1] "Total unique JMJD6 substrate proteins and sites:"
    ##    curated_oxK_site unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate               49           153
    ## 2:           Others             4172         49508
    ## [1] "by the availability of stocihiometry data"
    ##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate                    TRUE               47           147
    ## 2:  JMJD6_substrate                   FALSE                2             6
    ## 3:           Others                    TRUE             4172         49508
    ## [1] "by the data source"
    ##     curated_oxK_site stoichiometry_available           data_source
    ##  1:  JMJD6_substrate                   FALSE                      
    ##  2:  JMJD6_substrate                    TRUE         HEK293_WT_JQ1
    ##  3:  JMJD6_substrate                    TRUE HeLa_JMJD6FLAG_FLAGIP
    ##  4:  JMJD6_substrate                    TRUE    HeLa_JMJD6KO_J6pep
    ##  5:  JMJD6_substrate                    TRUE      HeLa_JMJD6KO_JQ1
    ##  6:  JMJD6_substrate                    TRUE         HeLa_WT_J6pep
    ##  7:  JMJD6_substrate                    TRUE           HeLa_WT_JQ1
    ##  8:  JMJD6_substrate                    TRUE           MCF7_WT_JQ1
    ##  9:           Others                    TRUE         HEK293_WT_JQ1
    ## 10:           Others                    TRUE HeLa_JMJD6FLAG_FLAGIP
    ## 11:           Others                    TRUE    HeLa_JMJD6KO_J6pep
    ## 12:           Others                    TRUE      HeLa_JMJD6KO_JQ1
    ## 13:           Others                    TRUE         HeLa_WT_J6pep
    ## 14:           Others                    TRUE           HeLa_WT_JQ1
    ## 15:           Others                    TRUE           MCF7_WT_JQ1
    ##     unique_protein_N unique_site_N
    ##  1:                2             6
    ##  2:                4            30
    ##  3:               36           117
    ##  4:               38            75
    ##  5:               18            57
    ##  6:               39            77
    ##  7:               22            68
    ##  8:                6            33
    ##  9:              242          1861
    ## 10:             1142         10968
    ## 11:             3692         39093
    ## 12:              678          9920
    ## 13:             3224         30053
    ## 14:              688          9642
    ## 15:              269          2459

``` r
print("Hydroxylation sites identified by non-unique peptides; see the main text.")
```

    ## [1] "Hydroxylation sites identified by non-unique peptides; see the main text."

``` r
stoichiometry.dt[
    order(Accession, position, -oxK_ratio)
][, head(.SD, 1), by = list(Accession, position)][
    curated_oxK_site == "JMJD6_substrate" &
    (stoichiometry_available == FALSE | oxK_ratio == 0),
    .(Accession, position, curated_oxK_site, screen)
]
```

    ##             Accession position curated_oxK_site    screen
    ## 1:  Q14331|FRG1_HUMAN       27  JMJD6_substrate FLAGJMJD6
    ## 2:  Q14331|FRG1_HUMAN       29  JMJD6_substrate FLAGJMJD6
    ## 3:  Q14331|FRG1_HUMAN       30  JMJD6_substrate FLAGJMJD6
    ## 4: Q9UQ35|SRRM2_HUMAN      241  JMJD6_substrate FLAGJMJD6
    ## 5: Q9UQ35|SRRM2_HUMAN      243  JMJD6_substrate FLAGJMJD6
    ## 6: Q9UQ35|SRRM2_HUMAN      244  JMJD6_substrate FLAGJMJD6

# Analysis of the properties of hydroxylation sites

``` r
stoichiometry.dt[, `:=`(
    protein_sequence = BSgenome::getSeq(
                             x = all.protein.bs,
                             name = Accession
                         ) %>% as.character
)]

seq.search.range <- 10

stoichiometry.dt[, `:=`(
    sequence = substr(
        protein_sequence,
        start = pmax(0, position - seq.search.range),
        stop = position + seq.search.range
    ),
    protein_sequence = NULL
)]

stoichiometry.dt[, full_length_seq_flag := (
    nchar(sequence) == (seq.search.range * 2 + 1)
)]

left.cols <- c(
    "data_source", "Accession", "uniprot_id", "position", "sequence",
    "full_length_seq_flag",
    "curated_oxK_site",
    "oxK_ratio"
)

stoichiometry.dt <- stoichiometry.dt[, c(
    left.cols, colnames(stoichiometry.dt)[!(colnames(stoichiometry.dt) %in% left.cols)]
), with = FALSE]

## Additional data filteration by the number of feature, and the availability of full length sequence

stoichiometry.dt <- stoichiometry.dt[
    curated_oxK_site == "JMJD6_substrate" |
    (
        residue == "K" &
        full_length_seq_flag == TRUE &
        total_n_feature_K > 2
    )
]

temp <- reportOxKSiteStats(stoichiometry.dt)
```

    ## [1] "Total unique JMJD6 substrate proteins and sites:"
    ##    curated_oxK_site unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate               49           153
    ## 2:           Others             2700         26311
    ## [1] "by the availability of stocihiometry data"
    ##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate                    TRUE               47           147
    ## 2:  JMJD6_substrate                   FALSE                2             6
    ## 3:           Others                    TRUE             2700         26311
    ## [1] "by the data source"
    ##     curated_oxK_site stoichiometry_available           data_source
    ##  1:  JMJD6_substrate                   FALSE                      
    ##  2:  JMJD6_substrate                    TRUE         HEK293_WT_JQ1
    ##  3:  JMJD6_substrate                    TRUE HeLa_JMJD6FLAG_FLAGIP
    ##  4:  JMJD6_substrate                    TRUE    HeLa_JMJD6KO_J6pep
    ##  5:  JMJD6_substrate                    TRUE      HeLa_JMJD6KO_JQ1
    ##  6:  JMJD6_substrate                    TRUE         HeLa_WT_J6pep
    ##  7:  JMJD6_substrate                    TRUE           HeLa_WT_JQ1
    ##  8:  JMJD6_substrate                    TRUE           MCF7_WT_JQ1
    ##  9:           Others                    TRUE         HEK293_WT_JQ1
    ## 10:           Others                    TRUE HeLa_JMJD6FLAG_FLAGIP
    ## 11:           Others                    TRUE    HeLa_JMJD6KO_J6pep
    ## 12:           Others                    TRUE      HeLa_JMJD6KO_JQ1
    ## 13:           Others                    TRUE         HeLa_WT_J6pep
    ## 14:           Others                    TRUE           HeLa_WT_JQ1
    ## 15:           Others                    TRUE           MCF7_WT_JQ1
    ##     unique_protein_N unique_site_N
    ##  1:                2             6
    ##  2:                4            30
    ##  3:               36           117
    ##  4:               38            75
    ##  5:               18            57
    ##  6:               39            77
    ##  7:               22            68
    ##  8:                6            33
    ##  9:               76           658
    ## 10:              907          5644
    ## 11:             2230         20913
    ## 12:              406          4473
    ## 13:             1966         15948
    ## 14:              437          4513
    ## 15:              128           943

# Analysis of sequence feature around hydroxylation site

## Amino acid enrichment around the hydroxylation sites

``` r
print("the following sites were excluded from the analysis:")
```

    ## [1] "the following sites were excluded from the analysis:"

``` r
stoichiometry.dt[
    (data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep")) |
    (curated_oxK_site == "JMJD6_substrate")
][CenterResidue != "K"][
    order(
        data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep"),
        curated_oxK_site == "JMJD6_substrate",
        oxK_ratio,
        decreasing = TRUE
    )
][
    !duplicated(paste(Accession_position))
]
```

    ##      data_source          Accession uniprot_id position        sequence
    ## 1: HeLa_WT_J6pep Q9Y5B9|SP16H_HUMAN     Q9Y5B9     1043 RGSRHSSAPPKKKRK
    ## 2:   HeLa_WT_JQ1 Q9Y5B9|SP16H_HUMAN     Q9Y5B9     1044  GSRHSSAPPKKKRK
    ##    full_length_seq_flag curated_oxK_site oxK_ratio    screen residue IUPRED2
    ## 1:                FALSE  JMJD6_substrate 0.2550824 JQ1,J6pep       K  0.7916
    ## 2:                FALSE  JMJD6_substrate 0.1028174 JQ1,J6pep       K  0.8235
    ##    K_position K_ratio K_ratio_score WindowHydropathy windowCharge CenterResidue
    ## 1:          1     0.4           0.4               NA           NA              
    ## 2:          1     0.3           0.4               NA           NA              
    ##    Window total_area_K total_area_oxK total_n_feature_K total_n_feature_oxK
    ## 1:           323750718       82583100                53                  21
    ## 2:            37340000        3839200                11                   2
    ##          seq5 MW_within_1 MW_within_2 stoichiometry_available
    ## 1: SSAPPKKKRK       FALSE       FALSE                    TRUE
    ## 2:  SAPPKKKRK       FALSE       FALSE                    TRUE
    ##         Accession_position
    ## 1: Q9Y5B9|SP16H_HUMAN_1043
    ## 2: Q9Y5B9|SP16H_HUMAN_1044

``` r
non.duplicated.stoichiometry.dt <- stoichiometry.dt[
    (data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep")) |
    (curated_oxK_site == "JMJD6_substrate")
][CenterResidue == "K"][
    order(
        data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep"),
        curated_oxK_site == "JMJD6_substrate",
        oxK_ratio,
        decreasing = TRUE
    )
][
    !duplicated(paste(Accession_position))
]

temp <- reportOxKSiteStats(non.duplicated.stoichiometry.dt)
```

    ## [1] "Total unique JMJD6 substrate proteins and sites:"
    ##    curated_oxK_site unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate               48           151
    ## 2:           Others             1992         16952
    ## [1] "by the availability of stocihiometry data"
    ##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate                    TRUE               46           145
    ## 2:  JMJD6_substrate                   FALSE                2             6
    ## 3:           Others                    TRUE             1992         16952
    ## [1] "by the data source"
    ##    curated_oxK_site stoichiometry_available           data_source
    ## 1:  JMJD6_substrate                   FALSE                      
    ## 2:  JMJD6_substrate                    TRUE HeLa_JMJD6FLAG_FLAGIP
    ## 3:  JMJD6_substrate                    TRUE         HeLa_WT_J6pep
    ## 4:  JMJD6_substrate                    TRUE           HeLa_WT_JQ1
    ## 5:           Others                    TRUE         HeLa_WT_J6pep
    ## 6:           Others                    TRUE           HeLa_WT_JQ1
    ##    unique_protein_N unique_site_N
    ## 1:                2             6
    ## 2:               12            39
    ## 3:               34            54
    ## 4:               13            52
    ## 5:             1938         12612
    ## 6:              437          4340

``` r
## Sanity checks
non.duplicated.stoichiometry.dt[duplicated(Accession_position)]
```

    ## Empty data.table (0 rows and 27 cols): data_source,Accession,uniprot_id,position,sequence,full_length_seq_flag...

``` r
non.duplicated.stoichiometry.dt[, table(CenterResidue)]
```

    ## CenterResidue
    ##     K 
    ## 17103

``` r
## Back to analysis
non.duplicated.stoichiometry.dt[
  , (paste0("position_", gsub("-", "m", seq(-10, 10)))) :=
        str_split(sequence, pattern = "")[[1]] %>% as.list,
    by = seq_len(nrow(non.duplicated.stoichiometry.dt))
]

key.id.cols <- c("data_source", "Accession", "position", "curated_oxK_site", "oxK_ratio")

m.non.duplicated.stoichiometry.dt <- melt(
    non.duplicated.stoichiometry.dt,
    id.vars = c(key.id.cols),
    measure.vars = paste0("position_", gsub("-", "m", seq(-10, 10))),
    variable.name = "position2KOH",
    value.name = "amino_acid"
)

aa.count.per.pos.dt <- m.non.duplicated.stoichiometry.dt[
  , .N, by = list(curated_oxK_site, position2KOH, amino_acid)
]

aa.count.per.pos.dt[
  , total_AA_per_pos := sum(N),
    by = list(curated_oxK_site, position2KOH)
]

aa.count.per.pos.dt[
  , odds := (N / total_AA_per_pos) %>% {./(1 - .)}
]

odds.ratio.dt <- dcast(
    aa.count.per.pos.dt,
    position2KOH + amino_acid ~ curated_oxK_site,
    value.var = c("odds", "N", "total_AA_per_pos")
)

setnafill(
    odds.ratio.dt, fill = 0, cols = c("odds_JMJD6_substrate", "N_JMJD6_substrate")
)

odds.ratio.dt[, total_AA_per_pos_JMJD6_substrate := case_when(
                    is.na(total_AA_per_pos_JMJD6_substrate) ~
                        mean(
                            total_AA_per_pos_JMJD6_substrate,
                            na.rm = TRUE
                        ) %>% as.integer,
                    TRUE ~ total_AA_per_pos_JMJD6_substrate
                )]

odds.ratio.dt[
   , binom_p := binom.test(
        x = N_JMJD6_substrate, n = total_AA_per_pos_JMJD6_substrate,
        p = N_Others / total_AA_per_pos_Others
    )$p.value,
    by = seq_len(nrow(odds.ratio.dt))
]


odds.ratio.dt[, `:=`(
    adj_binom_p = p.adjust(binom_p, method = "fdr"),
    log2_odds_ratio = log2(odds_JMJD6_substrate / odds_Others)
)]

odds.ratio.dt[
  , capped_log2_odds_ratio := case_when(
        is.infinite(log2_odds_ratio) ~ sign(log2_odds_ratio) * 5,
        TRUE ~ log2_odds_ratio
)]

aa.order.dt <-
    odds.ratio.dt[
      , list(sum_odds = sum(capped_log2_odds_ratio, na.rm = TRUE)), by = amino_acid
    ][
        order(sum_odds)
    ]

odds.ratio.dt[, `:=`(
    amino_acid = factor(amino_acid, levels = aa.order.dt[, amino_acid]),
    position = gsub("position_", "", position2KOH) %>%
        {gsub("m", "-", .)} %>%
        factor(levels = -10:10)
)]

ggplot(
    data = odds.ratio.dt,
    aes(
        x = position,
        y = amino_acid,
        fill = capped_log2_odds_ratio#,
        ##color = adj_binom_p < 0.1
    )
) +
    geom_tile(size = 0.85) +
    scale_fill_gradient2(name = "Odds ratio") +
    scale_color_manual(values = c("TRUE" = "purple", "FALSE" = NA), name = "p(adj) < 0.1") +
    ylab("Amino acid") + 
    theme(
        legend.position = "bottom"
    )
```

![](j3-properties-of-hydroxylation-sites_files/figure-gfm/amino_acid_enrichment-1.png)<!-- -->

## Biophysical properties

``` r
non.duplicated.stoichiometry.dt <-
    non.duplicated.stoichiometry.dt[
        data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep") &
        total_n_feature_K > 2 &
        !is.na(windowCharge)
    ]

temp <- reportOxKSiteStats(non.duplicated.stoichiometry.dt)
```

    ## [1] "Total unique JMJD6 substrate proteins and sites:"
    ##    curated_oxK_site unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate               35            84
    ## 2:           Others             1992         16952
    ## [1] "by the availability of stocihiometry data"
    ##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
    ## 1:  JMJD6_substrate                    TRUE               35            84
    ## 2:           Others                    TRUE             1992         16952
    ## [1] "by the data source"
    ##    curated_oxK_site stoichiometry_available   data_source unique_protein_N
    ## 1:  JMJD6_substrate                    TRUE HeLa_WT_J6pep               30
    ## 2:  JMJD6_substrate                    TRUE   HeLa_WT_JQ1                7
    ## 3:           Others                    TRUE HeLa_WT_J6pep             1938
    ## 4:           Others                    TRUE   HeLa_WT_JQ1              437
    ##    unique_site_N
    ## 1:            45
    ## 2:            39
    ## 3:         12612
    ## 4:          4340

``` r
## Sanity checks
non.duplicated.stoichiometry.dt[, table(CenterResidue)]
```

    ## CenterResidue
    ##     K 
    ## 17036

``` r
data.cols <- c("IUPRED2", "WindowHydropathy", "windowCharge", "K_ratio")

## non.duplicated.stoichiometry.dt[
##     curated_oxK_site == "JMJD6_substrate"
## ][order(windowCharge)][, .(
##      Accession, position, Window, windowCharge, IUPRED2
##  )]


ggplot(
    data = non.duplicated.stoichiometry.dt[order(curated_oxK_site != "Others")],
    aes(
        x = windowCharge,
        y = IUPRED2,
        size = ifelse(
            curated_oxK_site == "JMJD6_substrate",
            oxK_ratio + 0.01,
            0.01
        ),
        fill = curated_oxK_site,
        alpha = curated_oxK_site,
        shape = curated_oxK_site
    )
) +
    geom_jitter(width = 0.04, height = 0) +
    theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "vertical"
    ) +
    scale_size_continuous(
        name = "Hydroxylation rate [%]"
    ) +
    scale_alpha_manual(
        values = c(
            "JMJD6_substrate" = 1,
            "Others" = 0.075
        ),
        name = "JMJD6 substrate"
    ) +
    scale_shape_manual(
        values = c(
            "JMJD6_substrate" = "circle filled",
            "Others" = "circle small"
        ),
        name = "JMJD6 substrate"
    ) +
    scale_fill_manual(values = c(
                          "JMJD6_substrate" = "#EE6677",
                          "Others" = "gray60"
                      ), name = "JMJD6 substrate") +
    scale_color_manual(values = c(
                           "TRUE" = "black",
                           "FALSE" = "gray60"
                       ), name = "BRDs") +
    xlab("Local charge") +
    ylab("Disorderedness") +
    ggtitle("Biophysical propety around K")
```

![](j3-properties-of-hydroxylation-sites_files/figure-gfm/analysis_of_sequence_feature-1.png)<!-- -->

``` r
print("Test for disorderedness")
```

    ## [1] "Test for disorderedness"

``` r
non.duplicated.stoichiometry.dt %$%
    {wilcox.test(
         .[curated_oxK_site == "JMJD6_substrate", IUPRED2],
         .[curated_oxK_site == "Others", IUPRED2]
    )}
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  .[curated_oxK_site == "JMJD6_substrate", IUPRED2] and .[curated_oxK_site == "Others", IUPRED2]
    ## W = 1246636, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r
print("Test for local charge")
```

    ## [1] "Test for local charge"

``` r
non.duplicated.stoichiometry.dt %$%
    {wilcox.test(
         .[curated_oxK_site == "JMJD6_substrate", windowCharge],
         .[curated_oxK_site == "Others", windowCharge]
    )}
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  .[curated_oxK_site == "JMJD6_substrate", windowCharge] and .[curated_oxK_site == "Others", windowCharge]
    ## W = 1179628, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

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
    ##  package              * version  date       lib source        
    ##  Biobase                2.48.0   2020-04-27 [1] Bioconductor  
    ##  BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor  
    ##  BiocParallel           1.22.0   2020-04-27 [1] Bioconductor  
    ##  Biostrings           * 2.56.0   2020-04-27 [1] Bioconductor  
    ##  bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
    ##  BSgenome             * 1.56.0   2020-04-27 [1] Bioconductor  
    ##  cellranger             1.1.0    2016-07-27 [1] CRAN (R 4.0.0)
    ##  cli                    3.0.1    2021-07-17 [1] CRAN (R 4.0.0)
    ##  colorspace             1.4-1    2019-03-18 [1] CRAN (R 4.0.0)
    ##  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
    ##  data.table           * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)
    ##  DelayedArray           0.14.0   2020-04-27 [1] Bioconductor  
    ##  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
    ##  dplyr                * 1.0.0    2020-05-29 [1] CRAN (R 4.0.0)
    ##  ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
    ##  evaluate               0.14     2019-05-28 [1] CRAN (R 4.0.0)
    ##  farver                 2.0.3    2020-01-16 [1] CRAN (R 4.0.0)
    ##  generics               0.0.2    2018-11-29 [1] CRAN (R 4.0.0)
    ##  GenomeInfoDb         * 1.24.0   2020-04-27 [1] Bioconductor  
    ##  GenomeInfoDbData       1.2.3    2021-06-10 [1] Bioconductor  
    ##  GenomicAlignments      1.24.0   2020-04-27 [1] Bioconductor  
    ##  GenomicRanges        * 1.40.0   2020-04-27 [1] Bioconductor  
    ##  ggplot2              * 3.3.1    2020-05-28 [1] CRAN (R 4.0.0)
    ##  glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
    ##  gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.0.0)
    ##  htmltools              0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
    ##  IRanges              * 2.22.1   2020-04-28 [1] Bioconductor  
    ##  khroma               * 1.3.0    2019-10-26 [1] CRAN (R 4.0.0)
    ##  knitr                * 1.28     2020-02-06 [1] CRAN (R 4.0.0)
    ##  labeling               0.3      2014-08-23 [1] CRAN (R 4.0.0)
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
    ##  readxl               * 1.3.1    2019-03-13 [1] CRAN (R 4.0.0)
    ##  rlang                  0.4.11   2021-04-30 [1] CRAN (R 4.0.0)
    ##  rmarkdown            * 2.2      2020-05-31 [1] CRAN (R 4.0.0)
    ##  Rsamtools              2.4.0    2020-04-27 [1] Bioconductor  
    ##  rtracklayer          * 1.48.0   2020-04-27 [1] Bioconductor  
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
