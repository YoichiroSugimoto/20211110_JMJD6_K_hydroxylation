---
title: "Analysis of lysine hydroxylation stoichiometry and protein features"
author: "Yoichiro Sugimoto"
date: "16 December, 2021"
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


# Package import



```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)

library("Biostrings")
library("GenomicRanges")
library("BSgenome")
library("readxl")

processors <- 8

set.seed(1)
```



```r
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



```r
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


```r
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

ggplot(
    data = j6pep.substrate.dt,
    aes(
        x = accession_position,
        y = oxK_ratio * 100
    )
) +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim = c(0, 100)) +
    theme(
        axis.title.x = element_blank()
    ) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    ylab("Stoichiometry [%]")
```

![](j3-properties-of-hydroxylation-sites_files/figure-html/stoichiometry data-1.png)<!-- -->

```r
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

```
## [1] "Total unique JMJD6 substrate proteins and sites:"
##    curated_oxK_site unique_protein_N unique_site_N
## 1:  JMJD6_substrate               54           163
## 2:           Others             4172         49498
## [1] "by the availability of stocihiometry data"
##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
## 1:  JMJD6_substrate                    TRUE               52           157
## 2:  JMJD6_substrate                   FALSE                2             6
## 3:           Others                    TRUE             4172         49498
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
##  2:                4            31
##  3:               39           123
##  4:               41            82
##  5:               22            65
##  6:               42            84
##  7:               26            76
##  8:                7            35
##  9:              242          1860
## 10:             1142         10962
## 11:             3692         39086
## 12:              678          9912
## 13:             3224         30046
## 14:              688          9634
## 15:              269          2457
```

```r
stoichiometry.dt[
    curated_oxK_site == "JMJD6_substrate" & stoichiometry_available == FALSE,
    .(Accession, position, curated_oxK_site, screen)
]
```

```
##             Accession position curated_oxK_site    screen
## 1:  Q14331|FRG1_HUMAN       27  JMJD6_substrate FLAGJMJD6
## 2:  Q14331|FRG1_HUMAN       29  JMJD6_substrate FLAGJMJD6
## 3:  Q14331|FRG1_HUMAN       30  JMJD6_substrate FLAGJMJD6
## 4: Q9UQ35|SRRM2_HUMAN      241  JMJD6_substrate FLAGJMJD6
## 5: Q9UQ35|SRRM2_HUMAN      243  JMJD6_substrate FLAGJMJD6
## 6: Q9UQ35|SRRM2_HUMAN      244  JMJD6_substrate FLAGJMJD6
```

# Analysis of the properties of hydroxylation sites



```r
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

```
## [1] "Total unique JMJD6 substrate proteins and sites:"
##    curated_oxK_site unique_protein_N unique_site_N
## 1:  JMJD6_substrate               54           163
## 2:           Others             2700         26301
## [1] "by the availability of stocihiometry data"
##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
## 1:  JMJD6_substrate                    TRUE               52           157
## 2:  JMJD6_substrate                   FALSE                2             6
## 3:           Others                    TRUE             2700         26301
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
##  2:                4            31
##  3:               39           123
##  4:               41            82
##  5:               22            65
##  6:               42            84
##  7:               26            76
##  8:                7            35
##  9:               76           657
## 10:              907          5639
## 11:             2230         20906
## 12:              406          4465
## 13:             1966         15941
## 14:              437          4505
## 15:              128           942
```

# Analysis of sequence feature around hydroxylation site

## Amino acid enrichment around the hydroxylation sites



```r
print("the following sites were excluded from the analysis:")
```

```
## [1] "the following sites were excluded from the analysis:"
```

```r
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
```

```r
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

```
## [1] "Total unique JMJD6 substrate proteins and sites:"
##    curated_oxK_site unique_protein_N unique_site_N
## 1:  JMJD6_substrate               53           161
## 2:           Others             1992         16943
## [1] "by the availability of stocihiometry data"
##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
## 1:  JMJD6_substrate                    TRUE               51           155
## 2:  JMJD6_substrate                   FALSE                2             6
## 3:           Others                    TRUE             1992         16943
## [1] "by the data source"
##    curated_oxK_site stoichiometry_available           data_source
## 1:  JMJD6_substrate                   FALSE                      
## 2:  JMJD6_substrate                    TRUE HeLa_JMJD6FLAG_FLAGIP
## 3:  JMJD6_substrate                    TRUE         HeLa_WT_J6pep
## 4:  JMJD6_substrate                    TRUE           HeLa_WT_JQ1
## 5:  JMJD6_substrate                    TRUE           MCF7_WT_JQ1
## 6:           Others                    TRUE         HeLa_WT_J6pep
## 7:           Others                    TRUE           HeLa_WT_JQ1
##    unique_protein_N unique_site_N
## 1:                2             6
## 2:               12            39
## 3:               36            57
## 4:               15            58
## 5:                1             1
## 6:             1938         12609
## 7:              437          4334
```

```r
## Sanity checks
non.duplicated.stoichiometry.dt[duplicated(Accession_position)]
```

```
## Empty data.table (0 rows and 27 cols): data_source,Accession,uniprot_id,position,sequence,full_length_seq_flag...
```

```r
non.duplicated.stoichiometry.dt[, table(CenterResidue)]
```

```
## CenterResidue
##     K 
## 17104
```

```r
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

![](j3-properties-of-hydroxylation-sites_files/figure-html/amino_acid_enrichment-1.png)<!-- -->


## Biophysical properties



```r
non.duplicated.stoichiometry.dt <-
    non.duplicated.stoichiometry.dt[
        data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep") &
        total_n_feature_K > 2 &
        !is.na(windowCharge)
    ]

temp <- reportOxKSiteStats(non.duplicated.stoichiometry.dt)
```

```
## [1] "Total unique JMJD6 substrate proteins and sites:"
##    curated_oxK_site unique_protein_N unique_site_N
## 1:  JMJD6_substrate               39            93
## 2:           Others             1992         16943
## [1] "by the availability of stocihiometry data"
##    curated_oxK_site stoichiometry_available unique_protein_N unique_site_N
## 1:  JMJD6_substrate                    TRUE               39            93
## 2:           Others                    TRUE             1992         16943
## [1] "by the data source"
##    curated_oxK_site stoichiometry_available   data_source unique_protein_N
## 1:  JMJD6_substrate                    TRUE HeLa_WT_J6pep               32
## 2:  JMJD6_substrate                    TRUE   HeLa_WT_JQ1                9
## 3:           Others                    TRUE HeLa_WT_J6pep             1938
## 4:           Others                    TRUE   HeLa_WT_JQ1              437
##    unique_site_N
## 1:            48
## 2:            45
## 3:         12609
## 4:          4334
```

```r
## Sanity checks
non.duplicated.stoichiometry.dt[, table(CenterResidue)]
```

```
## CenterResidue
##     K 
## 17036
```

```r
data.cols <- c("IUPRED2", "WindowHydropathy", "windowCharge", "K_ratio")

non.duplicated.stoichiometry.dt[
    curated_oxK_site == "JMJD6_substrate"
][order(windowCharge)][, .(
     Accession, position, Window, windowCharge, IUPRED2
 )]
```

```
##              Accession position      Window  windowCharge IUPRED2
##  1: P26368|U2AF2_HUMAN       15 QLNENKQERDK -1.064196e-03  0.8162
##  2:  P02545|LMNA_HUMAN      341 RLLAEKEREMA -9.960256e-05  0.6712
##  3: P26368|U2AF2_HUMAN       70 LTRGAKEEHGG  3.444825e-03  0.8375
##  4: P18077|RL35A_HUMAN       45 EFYLGKRCAYV  7.050756e-02  0.1671
##  5:  O60885|BRD4_HUMAN      561 EVEENKKSKAK  8.841388e-02  0.9126
##  6:  O60885|BRD4_HUMAN      562 VEENKKSKAKE  8.841388e-02  0.8920
##  7:  O60885|BRD4_HUMAN      332 PVKPPKKDVPD  8.864616e-02  0.9433
##  8: O95232|LC7L3_HUMAN      388 GSDDKKSSVKS  8.864616e-02  0.8530
##  9: Q05519|SRS11_HUMAN      411 KDRERKSESDK  8.911182e-02  0.7916
## 10: Q05519|SRS11_HUMAN      242 IEPDKKEEKRR  8.927831e-02  0.7951
## 11:  Q8IWS0|PHF6_HUMAN      173 RKTNFKGLSED  8.961075e-02  0.7209
## 12:  P11387|TOP1_HUMAN       40 HKEHKKEKDRE  9.603473e-02  0.9287
## 13:  P11387|TOP1_HUMAN       25 DSHKHKDKHKD  9.924617e-02  0.9081
## 14:  P11387|TOP1_HUMAN       23 LNDSHKHKDKH  9.997793e-02  0.8872
## 15: O43143|DHX15_HUMAN       18 YPSGKKRAGTD  1.802729e-01  0.8530
## 16:  Q15059|BRD3_HUMAN      364 STVKRKMDGRE  1.805185e-01  0.7120
## 17:  Q13428|TCOF_HUMAN     1348 ESRKRKLSGDQ  1.805185e-01  0.8375
## 18: P0DMV8|HS71A_HUMAN      248 EEFKRKHKKDI  1.829324e-01  0.3774
## 19:  P11387|TOP1_HUMAN       36 REHRHKEHKKE  1.914515e-01  0.9190
## 20:  O60885|BRD4_HUMAN      544 KEKDKKEKKKE  2.679015e-01  0.8872
## 21:  P62851|RS25_HUMAN       14 KKDAGKSAKKD  2.688655e-01  0.7982
## 22:  O60885|BRD4_HUMAN      537 QQNKPKKKEKD  2.690320e-01  0.9106
## 23:  Q6PD62|CTR9_HUMAN      993 KAEKKKAPKPE  2.691985e-01  0.9649
## 24:  O60885|BRD4_HUMAN      535 QPQQNKPKKKE  2.697638e-01  0.8920
## 25:  P35251|RFC1_HUMAN       38 TLKAKKGIKEI  2.697638e-01  0.5707
## 26:  O60885|BRD4_HUMAN      572 EPPPKKTKKNN  2.697638e-01  0.9329
## 27:  Q15059|BRD3_HUMAN      651 SASGKKQAAKS  2.703290e-01  0.7250
## 28:  O60885|BRD4_HUMAN      291 KGVKRKADTTT  2.703954e-01  0.8623
## 29:  Q15059|BRD3_HUMAN      684 LSSSKKPARKE  2.705619e-01  0.8313
## 30: O95232|LC7L3_HUMAN      392 KKSSVKSGSRE  2.705619e-01  0.8681
## 31: Q13435|SF3B2_HUMAN      320 TVSVSKKEKNR  2.705619e-01  0.8198
## 32:  O60885|BRD4_HUMAN      552 KKEKHKRKEEV  2.732086e-01  0.8565
## 33: P11142|HSP7C_HUMAN      248 AEFKRKHKKDI  2.736074e-01  0.3321
## 34:  P83731|RL24_HUMAN       61 RRKHKKGQSEE  2.745720e-01  0.5456
## 35:  P11387|TOP1_HUMAN       39 RHKEHKKEKDR  2.776175e-01  0.9249
## 36:  O60885|BRD4_HUMAN      541 PKKKEKDKKEK  3.585764e-01  0.9287
## 37:  O60885|BRD4_HUMAN      538 QNKPKKKEKDK  3.591417e-01  0.9106
## 38:  P25440|BRD2_HUMAN      755 TKKPPKKANEK  3.598734e-01  0.8713
## 39:  O60885|BRD4_HUMAN      575 PKKTKKNNSSN  3.604387e-01  0.9168
## 40:  P25440|BRD2_HUMAN      589 PKKSKKASGSG  3.604387e-01  0.9562
## 41:  P25440|BRD2_HUMAN      586 PPQPKKSKKAS  3.604387e-01  0.9664
## 42: Q08945|SSRP1_HUMAN      524 QLKKAKMAKDR  3.605050e-01  0.6948
## 43:  Q14839|CHD4_HUMAN       67 DPKIPKSKRQK  3.605050e-01  0.8493
## 44:  O60885|BRD4_HUMAN      289 TKKGVKRKADT  3.605050e-01  0.9126
## 45:  Q2NL82|TSR1_HUMAN       42 LKTLSKKVRKE  3.606715e-01  0.5139
## 46:  Q8TDN6|BRX1_HUMAN       17 FAVQAKKPKRN  3.612368e-01  0.6227
## 47:  Q15059|BRD3_HUMAN      683 QLSSSKKPARK  3.612368e-01  0.8530
## 48:  Q15059|BRD3_HUMAN      643 QKKQRKPFSAS  3.612368e-01  0.6089
## 49: Q9NQ29|LUC7L_HUMAN      323 RDRSAKYKFSR  3.620884e-01  0.7718
## 50: Q9NQ29|LUC7L_HUMAN      325 RSAKYKFSRER  3.622549e-01  0.6851
## 51:  O60885|BRD4_HUMAN      547 DKKEKKKEKHK  3.623537e-01  0.8313
## 52:  O60885|BRD4_HUMAN      550 EKKKEKHKRKE  3.633183e-01  0.8681
## 53:  Q9BZ95|NSD3_HUMAN      199 KEKRKKSNKHD  3.637170e-01  0.8872
## 54:  O60885|BRD4_HUMAN      727 APKSKKKGHPG  3.642160e-01  0.9684
## 55: Q8N9Q2|SR1IP_HUMAN      142 KRKKEKHSSTP  3.644488e-01  0.8125
## 56:  Q9BVP2|GNL3_HUMAN       23 RYKIQKKVREH  3.652341e-01  0.7415
## 57: Q9Y383|LC7L2_HUMAN      266 SRSHSKNPKRS  3.658122e-01  0.9039
## 58:  Q9BVP2|GNL3_HUMAN       20 CHKRYKIQKKV  4.353016e-01  0.7951
## 59:  O60885|BRD4_HUMAN      539 NKPKKKEKDKK  4.492514e-01  0.9308
## 60:  Q15059|BRD3_HUMAN      489 PVNKPKKKKEK  4.499831e-01  0.8279
## 61:   P10412|H14_HUMAN      174 KAKSPKKAKAA  4.505484e-01  0.8596
## 62:   P10412|H14_HUMAN      180 KAKAAKPKKAP  4.505484e-01  0.9269
## 63:   P10412|H14_HUMAN      175 AKSPKKAKAAK  4.505484e-01  0.8920
## 64:  Q15059|BRD3_HUMAN      487 QAPVNKPKKKK  4.505484e-01  0.7718
## 65: Q96SB4|SRPK1_HUMAN       18 KRTKAKKDKAQ  4.506147e-01  0.7036
## 66:  Q9BVP2|GNL3_HUMAN      221 ITKRVKAKKNA  4.513465e-01  0.4441
## 67: Q96SB4|SRPK1_HUMAN       19 RTKAKKDKAQR  4.514128e-01  0.7951
## 68:  P25440|BRD2_HUMAN      546 PISKPKRKREK  4.515793e-01  0.8198
## 69:  Q12873|CHD3_HUMAN       55 RKRGPKKQKEN  4.515793e-01  0.8828
## 70: Q5BKY9|F133B_HUMAN       90 ESSSKKRQRKK  4.515793e-01  0.8013
## 71: Q8WXA9|SREK1_HUMAN      397 SPRTSKTIKRK  4.521446e-01  0.9211
## 72: Q8WXA9|SREK1_HUMAN      400 TSKTIKRKSSR  4.521446e-01  0.9329
## 73: Q8N9Q2|SR1IP_HUMAN      139 EKKKRKKEKHS  4.539932e-01  0.9190
## 74:  Q6UN15|FIP1_HUMAN      569 KHKKSKRSKEG  4.545585e-01  0.9230
## 75: Q9BRS2|RIOK1_HUMAN      555 IPKHVKKRKEK  4.545585e-01  0.7799
## 76: Q9Y383|LC7L2_HUMAN      269 HSKNPKRSRSR  4.567199e-01  0.9433
## 77: Q08945|SSRP1_HUMAN      521 KRKQLKKAKMA  5.414562e-01  0.7331
## 78:  O60885|BRD4_HUMAN      286 PVKTKKGVKRK  5.414562e-01  0.8894
## 79:   P62081|RS7_HUMAN      113 SRTKNKQKRPR  5.430524e-01  0.8655
## 80:  O60885|BRD4_HUMAN      548 KKEKKKEKHKR  5.441029e-01  0.8421
## 81: Q8N9Q2|SR1IP_HUMAN      140 KKKRKKEKHSS  5.446682e-01  0.9211
## 82: Q9BRS2|RIOK1_HUMAN      552 KNKIPKHVKKR  5.452334e-01  0.8162
## 83:  Q6UN15|FIP1_HUMAN      564 SHRRHKHKKSK  5.535860e-01  0.9013
## 84:   P62241|RS8_HUMAN      144 KKRSKKIQKKY  6.315530e-01  0.4652
## 85: Q66PJ3|AR6P4_HUMAN      294 GKYKDKRRKKK  6.316193e-01  0.8713
## 86: Q96SB4|SRPK1_HUMAN       16 RKKRTKAKKDK  6.316322e-01  0.5901
## 87: Q66PJ3|AR6P4_HUMAN      292 KRGKYKDKRRK  6.324174e-01  0.8681
## 88: Q66PJ3|AR6P4_HUMAN      290 RKKRGKYKDKR  6.324174e-01  0.8655
## 89:  Q6UN15|FIP1_HUMAN      567 RHKHKKSKRSK  6.399185e-01  0.9168
## 90:  Q6UN15|FIP1_HUMAN      566 RRHKHKKSKRS  6.407166e-01  0.8966
## 91:  P46100|ATRX_HUMAN     1424 YKQKKKRRRIK  7.232589e-01  0.6531
## 92: Q13435|SF3B2_HUMAN      332 RRNRKKKKKPQ  7.232717e-01  0.7672
## 93:  P46100|ATRX_HUMAN     1422 RSYKQKKKRRR  7.240570e-01  0.5992
##              Accession position      Window  windowCharge IUPRED2
```

```r
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
    ylab("Disordedness") +
    ggtitle("Biophysical propety around K")
```

![](j3-properties-of-hydroxylation-sites_files/figure-html/analysis of sequence feature-1.png)<!-- -->



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
```
