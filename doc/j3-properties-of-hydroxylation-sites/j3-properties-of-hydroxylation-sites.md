---
title: "Analysis of lysine hydroxylation stoichiometry and protein features"
author: "Yoichiro Sugimoto"
date: "08 December, 2021"
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


## Manually curated JMJD6 hydroxylation sites



```r
manually.curated.hydroxylation.dt <-
    file.path(
        j3.input.dir,
        "20210721_JMJD6_manually_curated_hydroxylation_site.xlsx"
    ) %>%
    read_excel %>%
    data.table

manually.curated.hydroxylation.dt <-
    manually.curated.hydroxylation.dt[, .(Accession, position)]

manually.curated.hydroxylation.dt[, `:=`(
    residue_at_curated_site =
        BSgenome::getSeq(
                      x = all.protein.bs,
                      name = Accession
                  ) %>% as.character %>%
        {substr(
            .,
            start = position,
            stop = position
        )},
    JMJD6_substrate_flag = "JMJD6_substrate"
)]

## Sanity check
manually.curated.hydroxylation.dt[, table(residue_at_curated_site)]
```

```
## residue_at_curated_site
##   K 
## 163
```

```r
print("Positions and proteins in manually curated data:")
```

```
## [1] "Positions and proteins in manually curated data:"
```

```r
manually.curated.hydroxylation.dt[, list(
    unique_protein_N = uniqueN(Accession),
    unique_site_N = uniqueN(paste0(Accession, "_", position))
)]
```

```
##    unique_protein_N unique_site_N
## 1:               54           163
```

## Stoichiometry data



```r
## cells <- c(
##     "HeLa_WT", "HeLa_J6KO", "HEK293", "MCF7", "HeLa_WT_J6PD", "HeLa_J6KO_J6PD"
## )

stoichiometry.dt <- fread(
    file.path(
        j2.1.res.dir,
        "long_K_stoichiometry_data.csv"
    )
)

stoichiometry.dt[, stoichiometry_available := TRUE]

stoichiometry.dt <- merge(
    manually.curated.hydroxylation.dt,
    stoichiometry.dt,
    all = TRUE,
    by = c("Accession", "position")
)

stoichiometry.dt[, `:=`(
    stoichiometry_available = case_when(
        is.na(stoichiometry_available) ~ FALSE,
        TRUE ~ TRUE
    ),
    JMJD6_substrate_flag = case_when(
        is.na(JMJD6_substrate_flag) ~ "Others",
        TRUE ~ JMJD6_substrate_flag
    ),
    Accession_position = paste(Accession, position, sep = "_")
)]

reportOxKSiteStats <- function(stoichiometry.dt){
    print("Positions and proteins in aggregated data:")
    stoichiometry.dt[, list(
        unique_protein_N = uniqueN(Accession),
        unique_site_N = uniqueN(Accession_position)
    ),
    by = list(
        JMJD6_substrate_flag,
        stoichiometry_available
    )
    ][order(JMJD6_substrate_flag, -stoichiometry_available)] %>%
        print

    stoichiometry.dt[, list(
        unique_protein_N = uniqueN(Accession),
        unique_site_N = uniqueN(Accession_position)
    ),
    by = list(
        JMJD6_substrate_flag,
        stoichiometry_available,
        cell
    )
    ][order(JMJD6_substrate_flag, cell)] %>%
        print

    return()
}

temp <- reportOxKSiteStats(stoichiometry.dt)
```

```
## [1] "Positions and proteins in aggregated data:"
##    JMJD6_substrate_flag stoichiometry_available unique_protein_N unique_site_N
## 1:      JMJD6_substrate                    TRUE               46           122
## 2:      JMJD6_substrate                   FALSE               13            41
## 3:               Others                    TRUE             4044         46332
##     JMJD6_substrate_flag stoichiometry_available           cell
##  1:      JMJD6_substrate                    TRUE         HEK293
##  2:      JMJD6_substrate                    TRUE      HeLa_J6KO
##  3:      JMJD6_substrate                    TRUE HeLa_J6KO_J6PD
##  4:      JMJD6_substrate                    TRUE        HeLa_WT
##  5:      JMJD6_substrate                    TRUE   HeLa_WT_J6PD
##  6:      JMJD6_substrate                    TRUE           MCF7
##  7:      JMJD6_substrate                   FALSE           <NA>
##  8:               Others                    TRUE         HEK293
##  9:               Others                    TRUE      HeLa_J6KO
## 10:               Others                    TRUE HeLa_J6KO_J6PD
## 11:               Others                    TRUE        HeLa_WT
## 12:               Others                    TRUE   HeLa_WT_J6PD
## 13:               Others                    TRUE           MCF7
##     unique_protein_N unique_site_N
##  1:                4            31
##  2:               22            65
##  3:               41            82
##  4:               26            76
##  5:               42            84
##  6:                6            34
##  7:               13            41
##  8:              242          1860
##  9:              678          9912
## 10:             3692         39086
## 11:              688          9634
## 12:             3224         30046
## 13:              269          2458
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
    "cell", "Accession", "uniprot_id", "position", "sequence",
    "full_length_seq_flag",
    "JMJD6_substrate_flag",
    "oxK_ratio"
)

stoichiometry.dt <- stoichiometry.dt[, c(
    left.cols, colnames(stoichiometry.dt)[!(colnames(stoichiometry.dt) %in% left.cols)]
), with = FALSE]

## Additional data filteration by the number of feature, and the availability of full length sequence

stoichiometry.dt <- stoichiometry.dt[
    JMJD6_substrate_flag == "JMJD6_substrate" |
    (
        residue == "K" &
        full_length_seq_flag == TRUE &
        total_n_feature_K > 2
    )
]

temp <- reportOxKSiteStats(stoichiometry.dt)
```

```
## [1] "Positions and proteins in aggregated data:"
##    JMJD6_substrate_flag stoichiometry_available unique_protein_N unique_site_N
## 1:      JMJD6_substrate                    TRUE               46           122
## 2:      JMJD6_substrate                   FALSE               13            41
## 3:               Others                    TRUE             2486         24249
##     JMJD6_substrate_flag stoichiometry_available           cell
##  1:      JMJD6_substrate                    TRUE         HEK293
##  2:      JMJD6_substrate                    TRUE      HeLa_J6KO
##  3:      JMJD6_substrate                    TRUE HeLa_J6KO_J6PD
##  4:      JMJD6_substrate                    TRUE        HeLa_WT
##  5:      JMJD6_substrate                    TRUE   HeLa_WT_J6PD
##  6:      JMJD6_substrate                    TRUE           MCF7
##  7:      JMJD6_substrate                   FALSE           <NA>
##  8:               Others                    TRUE         HEK293
##  9:               Others                    TRUE      HeLa_J6KO
## 10:               Others                    TRUE HeLa_J6KO_J6PD
## 11:               Others                    TRUE        HeLa_WT
## 12:               Others                    TRUE   HeLa_WT_J6PD
## 13:               Others                    TRUE           MCF7
##     unique_protein_N unique_site_N
##  1:                4            31
##  2:               22            65
##  3:               41            82
##  4:               26            76
##  5:               42            84
##  6:                6            34
##  7:               13            41
##  8:               76           657
##  9:              406          4465
## 10:             2230         20906
## 11:              437          4505
## 12:             1966         15941
## 13:              128           942
```

# Analysis of sequence feature around hydroxylation site

## Amino acid enrichment around the hydroxylation sites



```r
non.duplicated.stoichiometry.dt <- stoichiometry.dt[
    !(cell %in% c("HeLa_J6KO", "HeLa_J6KO_J6PD")) | is.na(cell)
][
    order(
        JMJD6_substrate_flag == "JMJD6_substrate",
        oxK_ratio,
        decreasing = TRUE
    )
][
    !duplicated(paste(Accession, position))
]

temp <- reportOxKSiteStats(non.duplicated.stoichiometry.dt)
```

```
## [1] "Positions and proteins in aggregated data:"
##    JMJD6_substrate_flag stoichiometry_available unique_protein_N unique_site_N
## 1:      JMJD6_substrate                    TRUE               46           117
## 2:      JMJD6_substrate                   FALSE               13            41
## 3:               Others                    TRUE             2012         17215
##    JMJD6_substrate_flag stoichiometry_available         cell unique_protein_N
## 1:      JMJD6_substrate                    TRUE       HEK293                3
## 2:      JMJD6_substrate                    TRUE      HeLa_WT               16
## 3:      JMJD6_substrate                    TRUE HeLa_WT_J6PD               37
## 4:      JMJD6_substrate                    TRUE         MCF7                4
## 5:      JMJD6_substrate                   FALSE         <NA>               13
## 6:               Others                    TRUE       HEK293               56
## 7:               Others                    TRUE      HeLa_WT              437
## 8:               Others                    TRUE HeLa_WT_J6PD             1928
## 9:               Others                    TRUE         MCF7               71
##    unique_site_N
## 1:            15
## 2:            36
## 3:            57
## 4:             9
## 5:            41
## 6:           252
## 7:          4313
## 8:         12463
## 9:           187
```

```r
## Sanity checks
print("To discuss #1")
```

```
## [1] "To discuss #1"
```

```r
non.duplicated.stoichiometry.dt[
    JMJD6_substrate_flag == "JMJD6_substrate" & oxK_ratio == 0
]
```

```
##             cell          Accession uniprot_id position              sequence
##  1: HeLa_WT_J6PD O15042|SR140_HUMAN     O15042      981 HKDSPRDVSKKAKRSPSGSRT
##  2:      HeLa_WT O95232|LC7L3_HUMAN     O95232      248 TEEPDRDERLKKEKQEREERE
##  3:      HeLa_WT O95232|LC7L3_HUMAN     O95232      249 EEPDRDERLKKEKQEREEREK
##  4:      HeLa_WT O95232|LC7L3_HUMAN     O95232      251 PDRDERLKKEKQEREEREKER
##  5: HeLa_WT_J6PD O95232|LC7L3_HUMAN     O95232      388 EKEKRGSDDKKSSVKSGSREK
##  6: HeLa_WT_J6PD O95232|LC7L3_HUMAN     O95232      392 RGSDDKKSSVKSGSREKQSED
##  7: HeLa_WT_J6PD  P02545|LMNA_HUMAN     P02545      341 RDTSRRLLAEKEREMAEMRAR
##  8:      HeLa_WT P0DMV8|HS71A_HUMAN     P0DMV8      248 VNHFVEEFKRKHKKDISQNKR
##  9:      HeLa_WT P11142|HSP7C_HUMAN     P11142      248 VNHFIAEFKRKHKKDISENKR
## 10:      HeLa_WT  P11387|TOP1_HUMAN     P11387       23 EADFRLNDSHKHKDKHKDREH
## 11:      HeLa_WT  P11387|TOP1_HUMAN     P11387       40 DREHRHKEHKKEKDREKSKHS
## 12:      HeLa_WT  P11387|TOP1_HUMAN     P11387      159 PKKIKTEDTKKEKKRKLEEEE
## 13: HeLa_WT_J6PD Q05519|SRS11_HUMAN     Q05519      411 SKDKEKDRERKSESDKDVKQV
## 14: HeLa_WT_J6PD  Q13428|TCOF_HUMAN     Q13428     1348 SRKGWESRKRKLSGDQPAART
## 15:         MCF7 Q13435|SF3B2_HUMAN     Q13435      320 ETEEDTVSVSKKEKNRKRRNR
## 16:      HeLa_WT Q13435|SF3B2_HUMAN     Q13435      331 KEKNRKRRNRKKKKKPQRVRG
## 17:      HeLa_WT Q13435|SF3B2_HUMAN     Q13435      334 NRKRRNRKKKKKPQRVRGVSS
## 18:      HeLa_WT Q13435|SF3B2_HUMAN     Q13435      335 RKRRNRKKKKKPQRVRGVSSE
## 19:      HeLa_WT Q14498|RBM39_HUMAN     Q14498      103 GRYRSPYSGPKFNSAIRGKIG
## 20:      HeLa_WT Q14498|RBM39_HUMAN     Q14498      119 RGKIGLPHSIKLSRRRSRSKS
## 21:      HeLa_WT  Q15059|BRD3_HUMAN     Q15059      643 VKSCLQKKQRKPFSASGKKQA
## 22: HeLa_WT_J6PD Q66PJ3|AR6P4_HUMAN     Q66PJ3      290 SSSDGRKKRGKYKDKRRKKKK
## 23: HeLa_WT_J6PD Q66PJ3|AR6P4_HUMAN     Q66PJ3      294 GRKKRGKYKDKRRKKKKKRKK
## 24:      HeLa_WT Q6NYC1|JMJD6_HUMAN     Q6NYC1      219 FPTSTPRELIKVTRDEGGNQQ
## 25: HeLa_WT_J6PD Q6NYC1|JMJD6_HUMAN     Q6NYC1      307 VWHKTVRGRPKLSRKWYRILK
## 26:      HeLa_WT  Q6UN15|FIP1_HUMAN     Q6UN15      564 SEEGDSHRRHKHKKSKRSKEG
## 27: HeLa_WT_J6PD  Q8IWS0|PHF6_HUMAN     Q8IWS0      173 RKGRPRKTNFKGLSEDTRSTS
## 28: HeLa_WT_J6PD Q8N9Q2|SR1IP_HUMAN     Q8N9Q2      142 KKEKKKRKKEKHSSTPNSSEF
## 29: HeLa_WT_J6PD Q8WXA9|SREK1_HUMAN     Q8WXA9      400 SRSPRTSKTIKRKSSRSPSPR
## 30: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       16 LALQARKKRTKAKKDKAQRKS
## 31: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       18 LQARKKRTKAKKDKAQRKSET
## 32: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       19 QARKKRTKAKKDKAQRKSETQ
## 33:      HeLa_WT Q9BRS2|RIOK1_HUMAN     Q9BRS2      535 TDPDIDKKERKKMVKEAQREK
## 34:      HeLa_WT Q9BRS2|RIOK1_HUMAN     Q9BRS2      539 IDKKERKKMVKEAQREKRKNK
## 35: HeLa_WT_J6PD Q9BRS2|RIOK1_HUMAN     Q9BRS2      555 KRKNKIPKHVKKRKEKTAKTK
## 36:      HeLa_WT  Q9BVP2|GNL3_HUMAN     Q9BVP2       20 SKRMTCHKRYKIQKKVREHHR
## 37:      HeLa_WT  Q9BVP2|GNL3_HUMAN     Q9BVP2       23 MTCHKRYKIQKKVREHHRKLR
## 38: HeLa_WT_J6PD Q9NQ29|LUC7L_HUMAN     Q9NQ29      323 HRRASRDRSAKYKFSRERASR
## 39: HeLa_WT_J6PD Q9NWB6|ARGL1_HUMAN     Q9NWB6       13 RSRSRSSSRSKHTKSSKHNKK
## 40: HeLa_WT_J6PD  Q9NYK5|RM39_HUMAN     Q9NYK5      322 IWDKLLERSRKMVTEDQSKAT
##             cell          Accession uniprot_id position              sequence
##     full_length_seq_flag JMJD6_substrate_flag oxK_ratio residue_at_curated_site
##  1:                 TRUE      JMJD6_substrate         0                       K
##  2:                 TRUE      JMJD6_substrate         0                       K
##  3:                 TRUE      JMJD6_substrate         0                       K
##  4:                 TRUE      JMJD6_substrate         0                       K
##  5:                 TRUE      JMJD6_substrate         0                       K
##  6:                 TRUE      JMJD6_substrate         0                       K
##  7:                 TRUE      JMJD6_substrate         0                       K
##  8:                 TRUE      JMJD6_substrate         0                       K
##  9:                 TRUE      JMJD6_substrate         0                       K
## 10:                 TRUE      JMJD6_substrate         0                       K
## 11:                 TRUE      JMJD6_substrate         0                       K
## 12:                 TRUE      JMJD6_substrate         0                       K
## 13:                 TRUE      JMJD6_substrate         0                       K
## 14:                 TRUE      JMJD6_substrate         0                       K
## 15:                 TRUE      JMJD6_substrate         0                       K
## 16:                 TRUE      JMJD6_substrate         0                       K
## 17:                 TRUE      JMJD6_substrate         0                       K
## 18:                 TRUE      JMJD6_substrate         0                       K
## 19:                 TRUE      JMJD6_substrate         0                       K
## 20:                 TRUE      JMJD6_substrate         0                       K
## 21:                 TRUE      JMJD6_substrate         0                       K
## 22:                 TRUE      JMJD6_substrate         0                       K
## 23:                 TRUE      JMJD6_substrate         0                       K
## 24:                 TRUE      JMJD6_substrate         0                       K
## 25:                 TRUE      JMJD6_substrate         0                       K
## 26:                 TRUE      JMJD6_substrate         0                       K
## 27:                 TRUE      JMJD6_substrate         0                       K
## 28:                 TRUE      JMJD6_substrate         0                       K
## 29:                 TRUE      JMJD6_substrate         0                       K
## 30:                 TRUE      JMJD6_substrate         0                       K
## 31:                 TRUE      JMJD6_substrate         0                       K
## 32:                 TRUE      JMJD6_substrate         0                       K
## 33:                 TRUE      JMJD6_substrate         0                       K
## 34:                 TRUE      JMJD6_substrate         0                       K
## 35:                 TRUE      JMJD6_substrate         0                       K
## 36:                 TRUE      JMJD6_substrate         0                       K
## 37:                 TRUE      JMJD6_substrate         0                       K
## 38:                 TRUE      JMJD6_substrate         0                       K
## 39:                 TRUE      JMJD6_substrate         0                       K
## 40:                 TRUE      JMJD6_substrate         0                       K
##     full_length_seq_flag JMJD6_substrate_flag oxK_ratio residue_at_curated_site
##     residue IUPRED2 K_position K_ratio K_ratio_score WindowHydropathy
##  1:       K  0.9606          1     0.2           0.3       0.31122222
##  2:       K  0.6851          1     0.3           0.3       0.17411111
##  3:       K  0.6806          1     0.2           0.3       0.17411111
##  4:       K  0.6755          1     0.2           0.3       0.17411111
##  5:       K  0.8530          1     0.2           0.4       0.29144444
##  6:       K  0.8681          1     0.2           0.3       0.35566667
##  7:       K  0.6712          1     0.1           0.1       0.40611111
##  8:       K  0.3774          1     0.4           0.4       0.16055556
##  9:       K  0.3321          1     0.4           0.4       0.16055556
## 10:       K  0.8872          1     0.4           0.4       0.13700000
## 11:       K  0.9287          1     0.4           0.5       0.08277778
## 12:       K  0.6755          1     0.4           0.5       0.11366667
## 13:       K  0.7916          1     0.3           0.4       0.14811111
## 14:       K  0.8375          1     0.1           0.3       0.27166667
## 15:       K  0.8198          1     0.4           0.4       0.35322222
## 16:       K  0.7982          1     0.5           0.6       0.04955556
## 17:       K  0.7799          1     0.2           0.6       0.08533333
## 18:       K  0.7718          1     0.1           0.6       0.19277778
## 19:       K  0.6576          1     0.2           0.2       0.40500000
## 20:       K  0.7629          1     0.2           0.2       0.36411111
## 21:       K  0.6089          1     0.3           0.4       0.28400000
## 22:       K  0.8655          1     0.6           0.6       0.13977778
## 23:       K  0.8713          1     0.7           0.7       0.08911111
## 24:       K  0.5707          1     0.1           0.1       0.40000000
## 25:       K  0.2399          1     0.2           0.2       0.24944444
## 26:       K  0.9013          1     0.5           0.5       0.11600000
## 27:       K  0.7209          1     0.1           0.3       0.37533333
## 28:       K  0.8125          1     0.1           0.7       0.18888889
## 29:       K  0.9329          1     0.2           0.3       0.31733333
## 30:       K  0.5901          1     0.5           0.6       0.17422222
## 31:       K  0.7036          1     0.4           0.6       0.24455556
## 32:       K  0.7951          1     0.3           0.6       0.25688889
## 33:       K  0.7799          1     0.3           0.5       0.23600000
## 34:       K  0.7459          1     0.3           0.5       0.31122222
## 35:       K  0.7799          1     0.5           0.5       0.20133333
## 36:       K  0.7951          1     0.3           0.4       0.20877778
## 37:       K  0.7415          1     0.3           0.4       0.30511111
## 38:       K  0.7718          1     0.2           0.2       0.32600000
## 39:       K  0.9649          1     0.4           0.4       0.26044444
## 40:       K  0.3356          1     0.2           0.2       0.31111111
##     residue IUPRED2 K_position K_ratio K_ratio_score WindowHydropathy
##     windowCharge total_area_K total_area_oxK total_n_feature_K
##  1:    0.3304832      3092200              0                 1
##  2:    0.1088334     15544000              0                 1
##  3:    0.1090369     15544000              0                 1
##  4:    0.1090369     15544000              0                 1
##  5:    0.1083453     35237900              0                19
##  6:    0.3313776     35237900              0                19
##  7:   -0.1112312     10320000              0                10
##  8:    0.3344090   1035400000              0                 5
##  9:    0.3344090   2206150200              0                 6
## 10:    0.1175786     24432980              0                15
## 11:    0.2235841      5188940              0                 5
## 12:    0.2189674       395630              0                 1
## 13:   -0.1113537      8528200              0                 7
## 14:    0.3314587      2533500              0                 8
## 15:    0.2195772       246829              0                 2
## 16:    0.8839988     78541000              0                 2
## 17:    0.7728893     78541000              0                 2
## 18:    0.6617798     78541000              0                 2
## 19:    0.1101184    104370000              0                 2
## 20:    0.3369697    185470000              0                 2
## 21:    0.4415117   4684755410              0                43
## 22:    0.5507356      8300500              0                 6
## 23:    0.6618451      8300500              0                 6
## 24:    0.1104997     10408000              0                 1
## 25:    0.5535966      2241800              0                 2
## 26:    0.5664711      7363900              0                 4
## 27:    0.1094432     91444930              0                17
## 28:    0.3353034     10117900              0                11
## 29:    0.4415117      7010810              0                 3
## 30:    0.5507513     36438930              0                13
## 31:    0.4406173     36829160              0                15
## 32:    0.3295078     36829160              0                15
## 33:    0.5509548     14011000              0                 1
## 34:    0.3306867     14011000              0                 1
## 35:    0.4454374     45521600              0                13
## 36:    0.5562467    121020000              0                 3
## 37:    0.3306710    121020000              0                 3
## 38:    0.2203335     14315270              0                11
## 39:    0.3359943       716380              0                 1
## 40:    0.1107032       441510              0                 1
##     windowCharge total_area_K total_area_oxK total_n_feature_K
##     total_n_feature_oxK        seq5 MW_within_1 MW_within_2
##  1:                   0 RDVSKKAKRSP       FALSE       FALSE
##  2:                   0 RDERLKKEKQE       FALSE       FALSE
##  3:                   0 DERLKKEKQER       FALSE       FALSE
##  4:                   0 RLKKEKQEREE       FALSE       FALSE
##  5:                   0 GSDDKKSSVKS       FALSE       FALSE
##  6:                   0 KKSSVKSGSRE       FALSE       FALSE
##  7:                   0 RLLAEKEREMA       FALSE       FALSE
##  8:                   0 EEFKRKHKKDI       FALSE       FALSE
##  9:                   0 AEFKRKHKKDI       FALSE       FALSE
## 10:                   0 LNDSHKHKDKH       FALSE       FALSE
## 11:                   0 HKEHKKEKDRE       FALSE       FALSE
## 12:                   0 TEDTKKEKKRK       FALSE       FALSE
## 13:                   0 KDRERKSESDK       FALSE       FALSE
## 14:                   0 ESRKRKLSGDQ       FALSE       FALSE
## 15:                   0 TVSVSKKEKNR       FALSE       FALSE
## 16:                   0 KRRNRKKKKKP       FALSE       FALSE
## 17:                   0 NRKKKKKPQRV       FALSE       FALSE
## 18:                   0 RKKKKKPQRVR       FALSE       FALSE
## 19:                   0 PYSGPKFNSAI       FALSE       FALSE
## 20:                   0 LPHSIKLSRRR       FALSE       FALSE
## 21:                   0 QKKQRKPFSAS       FALSE       FALSE
## 22:                   0 RKKRGKYKDKR       FALSE       FALSE
## 23:                   0 GKYKDKRRKKK       FALSE       FALSE
## 24:                   0 PRELIKVTRDE       FALSE       FALSE
## 25:                   0 VRGRPKLSRKW       FALSE       FALSE
## 26:                   0 SHRRHKHKKSK       FALSE       FALSE
## 27:                   0 RKTNFKGLSED       FALSE       FALSE
## 28:                   0 KRKKEKHSSTP       FALSE       FALSE
## 29:                   0 TSKTIKRKSSR       FALSE       FALSE
## 30:                   0 RKKRTKAKKDK       FALSE       FALSE
## 31:                   0 KRTKAKKDKAQ       FALSE       FALSE
## 32:                   0 RTKAKKDKAQR       FALSE       FALSE
## 33:                   0 DKKERKKMVKE       FALSE        TRUE
## 34:                   0 RKKMVKEAQRE       FALSE       FALSE
## 35:                   0 IPKHVKKRKEK       FALSE       FALSE
## 36:                   0 CHKRYKIQKKV       FALSE       FALSE
## 37:                   0 RYKIQKKVREH       FALSE       FALSE
## 38:                   0 RDRSAKYKFSR       FALSE       FALSE
## 39:                   0 SSSRSKHTKSS       FALSE       FALSE
## 40:                   0 LERSRKMVTED        TRUE        TRUE
##     total_n_feature_oxK        seq5 MW_within_1 MW_within_2
##     stoichiometry_available     Accession_position
##  1:                    TRUE O15042|SR140_HUMAN_981
##  2:                    TRUE O95232|LC7L3_HUMAN_248
##  3:                    TRUE O95232|LC7L3_HUMAN_249
##  4:                    TRUE O95232|LC7L3_HUMAN_251
##  5:                    TRUE O95232|LC7L3_HUMAN_388
##  6:                    TRUE O95232|LC7L3_HUMAN_392
##  7:                    TRUE  P02545|LMNA_HUMAN_341
##  8:                    TRUE P0DMV8|HS71A_HUMAN_248
##  9:                    TRUE P11142|HSP7C_HUMAN_248
## 10:                    TRUE   P11387|TOP1_HUMAN_23
## 11:                    TRUE   P11387|TOP1_HUMAN_40
## 12:                    TRUE  P11387|TOP1_HUMAN_159
## 13:                    TRUE Q05519|SRS11_HUMAN_411
## 14:                    TRUE Q13428|TCOF_HUMAN_1348
## 15:                    TRUE Q13435|SF3B2_HUMAN_320
## 16:                    TRUE Q13435|SF3B2_HUMAN_331
## 17:                    TRUE Q13435|SF3B2_HUMAN_334
## 18:                    TRUE Q13435|SF3B2_HUMAN_335
## 19:                    TRUE Q14498|RBM39_HUMAN_103
## 20:                    TRUE Q14498|RBM39_HUMAN_119
## 21:                    TRUE  Q15059|BRD3_HUMAN_643
## 22:                    TRUE Q66PJ3|AR6P4_HUMAN_290
## 23:                    TRUE Q66PJ3|AR6P4_HUMAN_294
## 24:                    TRUE Q6NYC1|JMJD6_HUMAN_219
## 25:                    TRUE Q6NYC1|JMJD6_HUMAN_307
## 26:                    TRUE  Q6UN15|FIP1_HUMAN_564
## 27:                    TRUE  Q8IWS0|PHF6_HUMAN_173
## 28:                    TRUE Q8N9Q2|SR1IP_HUMAN_142
## 29:                    TRUE Q8WXA9|SREK1_HUMAN_400
## 30:                    TRUE  Q96SB4|SRPK1_HUMAN_16
## 31:                    TRUE  Q96SB4|SRPK1_HUMAN_18
## 32:                    TRUE  Q96SB4|SRPK1_HUMAN_19
## 33:                    TRUE Q9BRS2|RIOK1_HUMAN_535
## 34:                    TRUE Q9BRS2|RIOK1_HUMAN_539
## 35:                    TRUE Q9BRS2|RIOK1_HUMAN_555
## 36:                    TRUE   Q9BVP2|GNL3_HUMAN_20
## 37:                    TRUE   Q9BVP2|GNL3_HUMAN_23
## 38:                    TRUE Q9NQ29|LUC7L_HUMAN_323
## 39:                    TRUE  Q9NWB6|ARGL1_HUMAN_13
## 40:                    TRUE  Q9NYK5|RM39_HUMAN_322
##     stoichiometry_available     Accession_position
```

```r
non.duplicated.stoichiometry.dt[duplicated(Accession_position)]
```

```
## Empty data.table (0 rows and 25 cols): cell,Accession,uniprot_id,position,sequence,full_length_seq_flag...
```

```r
## Back to analysis
non.duplicated.stoichiometry.dt[
  , (paste0("position_", gsub("-", "m", seq(-10, 10)))) :=
        str_split(sequence, pattern = "")[[1]] %>% as.list,
    by = seq_len(nrow(non.duplicated.stoichiometry.dt))
]

key.id.cols <- c("cell", "Accession", "position", "JMJD6_substrate_flag", "oxK_ratio")

m.non.duplicated.stoichiometry.dt <- melt(
    non.duplicated.stoichiometry.dt,
    id.vars = c(key.id.cols),
    measure.vars = paste0("position_", gsub("-", "m", seq(-10, 10))),
    variable.name = "position2KOH",
    value.name = "amino_acid"
)

aa.count.per.pos.dt <- m.non.duplicated.stoichiometry.dt[
  , .N, by = list(JMJD6_substrate_flag, position2KOH, amino_acid)
]

aa.count.per.pos.dt[
  , total_AA_per_pos := sum(N),
    by = list(JMJD6_substrate_flag, position2KOH)
]

aa.count.per.pos.dt[
  , odds := (N / total_AA_per_pos) %>% {./(1 - .)}
]

odds.ratio.dt <- dcast(
    aa.count.per.pos.dt,
    position2KOH + amino_acid ~ JMJD6_substrate_flag,
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

![](j3-properties-of-hydroxylation-sites_files/figure-html/amino acid enrichment-1.png)<!-- -->


## Biophysical properties



```r
non.duplicated.stoichiometry.dt <-
    non.duplicated.stoichiometry.dt[
        total_n_feature_K > 2
    ]

## Sanity checks
print("To discuss #1")
```

```
## [1] "To discuss #1"
```

```r
non.duplicated.stoichiometry.dt[
    JMJD6_substrate_flag == "JMJD6_substrate" & oxK_ratio == 0
]
```

```
##             cell          Accession uniprot_id position              sequence
##  1: HeLa_WT_J6PD O95232|LC7L3_HUMAN     O95232      388 EKEKRGSDDKKSSVKSGSREK
##  2: HeLa_WT_J6PD O95232|LC7L3_HUMAN     O95232      392 RGSDDKKSSVKSGSREKQSED
##  3: HeLa_WT_J6PD  P02545|LMNA_HUMAN     P02545      341 RDTSRRLLAEKEREMAEMRAR
##  4:      HeLa_WT P0DMV8|HS71A_HUMAN     P0DMV8      248 VNHFVEEFKRKHKKDISQNKR
##  5:      HeLa_WT P11142|HSP7C_HUMAN     P11142      248 VNHFIAEFKRKHKKDISENKR
##  6:      HeLa_WT  P11387|TOP1_HUMAN     P11387       23 EADFRLNDSHKHKDKHKDREH
##  7:      HeLa_WT  P11387|TOP1_HUMAN     P11387       40 DREHRHKEHKKEKDREKSKHS
##  8: HeLa_WT_J6PD Q05519|SRS11_HUMAN     Q05519      411 SKDKEKDRERKSESDKDVKQV
##  9: HeLa_WT_J6PD  Q13428|TCOF_HUMAN     Q13428     1348 SRKGWESRKRKLSGDQPAART
## 10:      HeLa_WT  Q15059|BRD3_HUMAN     Q15059      643 VKSCLQKKQRKPFSASGKKQA
## 11: HeLa_WT_J6PD Q66PJ3|AR6P4_HUMAN     Q66PJ3      290 SSSDGRKKRGKYKDKRRKKKK
## 12: HeLa_WT_J6PD Q66PJ3|AR6P4_HUMAN     Q66PJ3      294 GRKKRGKYKDKRRKKKKKRKK
## 13:      HeLa_WT  Q6UN15|FIP1_HUMAN     Q6UN15      564 SEEGDSHRRHKHKKSKRSKEG
## 14: HeLa_WT_J6PD  Q8IWS0|PHF6_HUMAN     Q8IWS0      173 RKGRPRKTNFKGLSEDTRSTS
## 15: HeLa_WT_J6PD Q8N9Q2|SR1IP_HUMAN     Q8N9Q2      142 KKEKKKRKKEKHSSTPNSSEF
## 16: HeLa_WT_J6PD Q8WXA9|SREK1_HUMAN     Q8WXA9      400 SRSPRTSKTIKRKSSRSPSPR
## 17: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       16 LALQARKKRTKAKKDKAQRKS
## 18: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       18 LQARKKRTKAKKDKAQRKSET
## 19: HeLa_WT_J6PD Q96SB4|SRPK1_HUMAN     Q96SB4       19 QARKKRTKAKKDKAQRKSETQ
## 20: HeLa_WT_J6PD Q9BRS2|RIOK1_HUMAN     Q9BRS2      555 KRKNKIPKHVKKRKEKTAKTK
## 21:      HeLa_WT  Q9BVP2|GNL3_HUMAN     Q9BVP2       20 SKRMTCHKRYKIQKKVREHHR
## 22:      HeLa_WT  Q9BVP2|GNL3_HUMAN     Q9BVP2       23 MTCHKRYKIQKKVREHHRKLR
## 23: HeLa_WT_J6PD Q9NQ29|LUC7L_HUMAN     Q9NQ29      323 HRRASRDRSAKYKFSRERASR
##             cell          Accession uniprot_id position              sequence
##     full_length_seq_flag JMJD6_substrate_flag oxK_ratio residue_at_curated_site
##  1:                 TRUE      JMJD6_substrate         0                       K
##  2:                 TRUE      JMJD6_substrate         0                       K
##  3:                 TRUE      JMJD6_substrate         0                       K
##  4:                 TRUE      JMJD6_substrate         0                       K
##  5:                 TRUE      JMJD6_substrate         0                       K
##  6:                 TRUE      JMJD6_substrate         0                       K
##  7:                 TRUE      JMJD6_substrate         0                       K
##  8:                 TRUE      JMJD6_substrate         0                       K
##  9:                 TRUE      JMJD6_substrate         0                       K
## 10:                 TRUE      JMJD6_substrate         0                       K
## 11:                 TRUE      JMJD6_substrate         0                       K
## 12:                 TRUE      JMJD6_substrate         0                       K
## 13:                 TRUE      JMJD6_substrate         0                       K
## 14:                 TRUE      JMJD6_substrate         0                       K
## 15:                 TRUE      JMJD6_substrate         0                       K
## 16:                 TRUE      JMJD6_substrate         0                       K
## 17:                 TRUE      JMJD6_substrate         0                       K
## 18:                 TRUE      JMJD6_substrate         0                       K
## 19:                 TRUE      JMJD6_substrate         0                       K
## 20:                 TRUE      JMJD6_substrate         0                       K
## 21:                 TRUE      JMJD6_substrate         0                       K
## 22:                 TRUE      JMJD6_substrate         0                       K
## 23:                 TRUE      JMJD6_substrate         0                       K
##     full_length_seq_flag JMJD6_substrate_flag oxK_ratio residue_at_curated_site
##     residue IUPRED2 K_position K_ratio K_ratio_score WindowHydropathy
##  1:       K  0.8530          1     0.2           0.4       0.29144444
##  2:       K  0.8681          1     0.2           0.3       0.35566667
##  3:       K  0.6712          1     0.1           0.1       0.40611111
##  4:       K  0.3774          1     0.4           0.4       0.16055556
##  5:       K  0.3321          1     0.4           0.4       0.16055556
##  6:       K  0.8872          1     0.4           0.4       0.13700000
##  7:       K  0.9287          1     0.4           0.5       0.08277778
##  8:       K  0.7916          1     0.3           0.4       0.14811111
##  9:       K  0.8375          1     0.1           0.3       0.27166667
## 10:       K  0.6089          1     0.3           0.4       0.28400000
## 11:       K  0.8655          1     0.6           0.6       0.13977778
## 12:       K  0.8713          1     0.7           0.7       0.08911111
## 13:       K  0.9013          1     0.5           0.5       0.11600000
## 14:       K  0.7209          1     0.1           0.3       0.37533333
## 15:       K  0.8125          1     0.1           0.7       0.18888889
## 16:       K  0.9329          1     0.2           0.3       0.31733333
## 17:       K  0.5901          1     0.5           0.6       0.17422222
## 18:       K  0.7036          1     0.4           0.6       0.24455556
## 19:       K  0.7951          1     0.3           0.6       0.25688889
## 20:       K  0.7799          1     0.5           0.5       0.20133333
## 21:       K  0.7951          1     0.3           0.4       0.20877778
## 22:       K  0.7415          1     0.3           0.4       0.30511111
## 23:       K  0.7718          1     0.2           0.2       0.32600000
##     residue IUPRED2 K_position K_ratio K_ratio_score WindowHydropathy
##     windowCharge total_area_K total_area_oxK total_n_feature_K
##  1:    0.1083453     35237900              0                19
##  2:    0.3313776     35237900              0                19
##  3:   -0.1112312     10320000              0                10
##  4:    0.3344090   1035400000              0                 5
##  5:    0.3344090   2206150200              0                 6
##  6:    0.1175786     24432980              0                15
##  7:    0.2235841      5188940              0                 5
##  8:   -0.1113537      8528200              0                 7
##  9:    0.3314587      2533500              0                 8
## 10:    0.4415117   4684755410              0                43
## 11:    0.5507356      8300500              0                 6
## 12:    0.6618451      8300500              0                 6
## 13:    0.5664711      7363900              0                 4
## 14:    0.1094432     91444930              0                17
## 15:    0.3353034     10117900              0                11
## 16:    0.4415117      7010810              0                 3
## 17:    0.5507513     36438930              0                13
## 18:    0.4406173     36829160              0                15
## 19:    0.3295078     36829160              0                15
## 20:    0.4454374     45521600              0                13
## 21:    0.5562467    121020000              0                 3
## 22:    0.3306710    121020000              0                 3
## 23:    0.2203335     14315270              0                11
##     windowCharge total_area_K total_area_oxK total_n_feature_K
##     total_n_feature_oxK        seq5 MW_within_1 MW_within_2
##  1:                   0 GSDDKKSSVKS       FALSE       FALSE
##  2:                   0 KKSSVKSGSRE       FALSE       FALSE
##  3:                   0 RLLAEKEREMA       FALSE       FALSE
##  4:                   0 EEFKRKHKKDI       FALSE       FALSE
##  5:                   0 AEFKRKHKKDI       FALSE       FALSE
##  6:                   0 LNDSHKHKDKH       FALSE       FALSE
##  7:                   0 HKEHKKEKDRE       FALSE       FALSE
##  8:                   0 KDRERKSESDK       FALSE       FALSE
##  9:                   0 ESRKRKLSGDQ       FALSE       FALSE
## 10:                   0 QKKQRKPFSAS       FALSE       FALSE
## 11:                   0 RKKRGKYKDKR       FALSE       FALSE
## 12:                   0 GKYKDKRRKKK       FALSE       FALSE
## 13:                   0 SHRRHKHKKSK       FALSE       FALSE
## 14:                   0 RKTNFKGLSED       FALSE       FALSE
## 15:                   0 KRKKEKHSSTP       FALSE       FALSE
## 16:                   0 TSKTIKRKSSR       FALSE       FALSE
## 17:                   0 RKKRTKAKKDK       FALSE       FALSE
## 18:                   0 KRTKAKKDKAQ       FALSE       FALSE
## 19:                   0 RTKAKKDKAQR       FALSE       FALSE
## 20:                   0 IPKHVKKRKEK       FALSE       FALSE
## 21:                   0 CHKRYKIQKKV       FALSE       FALSE
## 22:                   0 RYKIQKKVREH       FALSE       FALSE
## 23:                   0 RDRSAKYKFSR       FALSE       FALSE
##     total_n_feature_oxK        seq5 MW_within_1 MW_within_2
##     stoichiometry_available     Accession_position position_m10 position_m9
##  1:                    TRUE O95232|LC7L3_HUMAN_388            E           K
##  2:                    TRUE O95232|LC7L3_HUMAN_392            R           G
##  3:                    TRUE  P02545|LMNA_HUMAN_341            R           D
##  4:                    TRUE P0DMV8|HS71A_HUMAN_248            V           N
##  5:                    TRUE P11142|HSP7C_HUMAN_248            V           N
##  6:                    TRUE   P11387|TOP1_HUMAN_23            E           A
##  7:                    TRUE   P11387|TOP1_HUMAN_40            D           R
##  8:                    TRUE Q05519|SRS11_HUMAN_411            S           K
##  9:                    TRUE Q13428|TCOF_HUMAN_1348            S           R
## 10:                    TRUE  Q15059|BRD3_HUMAN_643            V           K
## 11:                    TRUE Q66PJ3|AR6P4_HUMAN_290            S           S
## 12:                    TRUE Q66PJ3|AR6P4_HUMAN_294            G           R
## 13:                    TRUE  Q6UN15|FIP1_HUMAN_564            S           E
## 14:                    TRUE  Q8IWS0|PHF6_HUMAN_173            R           K
## 15:                    TRUE Q8N9Q2|SR1IP_HUMAN_142            K           K
## 16:                    TRUE Q8WXA9|SREK1_HUMAN_400            S           R
## 17:                    TRUE  Q96SB4|SRPK1_HUMAN_16            L           A
## 18:                    TRUE  Q96SB4|SRPK1_HUMAN_18            L           Q
## 19:                    TRUE  Q96SB4|SRPK1_HUMAN_19            Q           A
## 20:                    TRUE Q9BRS2|RIOK1_HUMAN_555            K           R
## 21:                    TRUE   Q9BVP2|GNL3_HUMAN_20            S           K
## 22:                    TRUE   Q9BVP2|GNL3_HUMAN_23            M           T
## 23:                    TRUE Q9NQ29|LUC7L_HUMAN_323            H           R
##     stoichiometry_available     Accession_position position_m10 position_m9
##     position_m8 position_m7 position_m6 position_m5 position_m4 position_m3
##  1:           E           K           R           G           S           D
##  2:           S           D           D           K           K           S
##  3:           T           S           R           R           L           L
##  4:           H           F           V           E           E           F
##  5:           H           F           I           A           E           F
##  6:           D           F           R           L           N           D
##  7:           E           H           R           H           K           E
##  8:           D           K           E           K           D           R
##  9:           K           G           W           E           S           R
## 10:           S           C           L           Q           K           K
## 11:           S           D           G           R           K           K
## 12:           K           K           R           G           K           Y
## 13:           E           G           D           S           H           R
## 14:           G           R           P           R           K           T
## 15:           E           K           K           K           R           K
## 16:           S           P           R           T           S           K
## 17:           L           Q           A           R           K           K
## 18:           A           R           K           K           R           T
## 19:           R           K           K           R           T           K
## 20:           K           N           K           I           P           K
## 21:           R           M           T           C           H           K
## 22:           C           H           K           R           Y           K
## 23:           R           A           S           R           D           R
##     position_m8 position_m7 position_m6 position_m5 position_m4 position_m3
##     position_m2 position_m1 position_0 position_1 position_2 position_3
##  1:           D           K          K          S          S          V
##  2:           S           V          K          S          G          S
##  3:           A           E          K          E          R          E
##  4:           K           R          K          H          K          K
##  5:           K           R          K          H          K          K
##  6:           S           H          K          H          K          D
##  7:           H           K          K          E          K          D
##  8:           E           R          K          S          E          S
##  9:           K           R          K          L          S          G
## 10:           Q           R          K          P          F          S
## 11:           R           G          K          Y          K          D
## 12:           K           D          K          R          R          K
## 13:           R           H          K          H          K          K
## 14:           N           F          K          G          L          S
## 15:           K           E          K          H          S          S
## 16:           T           I          K          R          K          S
## 17:           R           T          K          A          K          K
## 18:           K           A          K          K          D          K
## 19:           A           K          K          D          K          A
## 20:           H           V          K          K          R          K
## 21:           R           Y          K          I          Q          K
## 22:           I           Q          K          K          V          R
## 23:           S           A          K          Y          K          F
##     position_m2 position_m1 position_0 position_1 position_2 position_3
##     position_4 position_5 position_6 position_7 position_8 position_9
##  1:          K          S          G          S          R          E
##  2:          R          E          K          Q          S          E
##  3:          M          A          E          M          R          A
##  4:          D          I          S          Q          N          K
##  5:          D          I          S          E          N          K
##  6:          K          H          K          D          R          E
##  7:          R          E          K          S          K          H
##  8:          D          K          D          V          K          Q
##  9:          D          Q          P          A          A          R
## 10:          A          S          G          K          K          Q
## 11:          K          R          R          K          K          K
## 12:          K          K          K          K          R          K
## 13:          S          K          R          S          K          E
## 14:          E          D          T          R          S          T
## 15:          T          P          N          S          S          E
## 16:          S          R          S          P          S          P
## 17:          D          K          A          Q          R          K
## 18:          A          Q          R          K          S          E
## 19:          Q          R          K          S          E          T
## 20:          E          K          T          A          K          T
## 21:          K          V          R          E          H          H
## 22:          E          H          H          R          K          L
## 23:          S          R          E          R          A          S
##     position_4 position_5 position_6 position_7 position_8 position_9
##     position_10
##  1:           K
##  2:           D
##  3:           R
##  4:           R
##  5:           R
##  6:           H
##  7:           S
##  8:           V
##  9:           T
## 10:           A
## 11:           K
## 12:           K
## 13:           G
## 14:           S
## 15:           F
## 16:           R
## 17:           S
## 18:           T
## 19:           Q
## 20:           K
## 21:           R
## 22:           R
## 23:           R
##     position_10
```

```r
data.cols <- c("IUPRED2", "WindowHydropathy", "windowCharge", "K_ratio")

m.stoichiometry.dt <- melt(
    non.duplicated.stoichiometry.dt,
    measure.vars = data.cols,
    id.vars = c("cell", "Accession", "position", "JMJD6_substrate_flag", "oxK_ratio")
)

temp <- reportOxKSiteStats(non.duplicated.stoichiometry.dt)
```

```
## [1] "Positions and proteins in aggregated data:"
##    JMJD6_substrate_flag stoichiometry_available unique_protein_N unique_site_N
## 1:      JMJD6_substrate                    TRUE               40            94
## 2:               Others                    TRUE             2012         17215
##    JMJD6_substrate_flag stoichiometry_available         cell unique_protein_N
## 1:      JMJD6_substrate                    TRUE       HEK293                3
## 2:      JMJD6_substrate                    TRUE      HeLa_WT               10
## 3:      JMJD6_substrate                    TRUE HeLa_WT_J6PD               33
## 4:      JMJD6_substrate                    TRUE         MCF7                3
## 5:               Others                    TRUE       HEK293               56
## 6:               Others                    TRUE      HeLa_WT              437
## 7:               Others                    TRUE HeLa_WT_J6PD             1928
## 8:               Others                    TRUE         MCF7               71
##    unique_site_N
## 1:            15
## 2:            23
## 3:            48
## 4:             8
## 5:           252
## 6:          4313
## 7:         12463
## 8:           187
```

```r
ggplot(
    data = non.duplicated.stoichiometry.dt[order(JMJD6_substrate_flag != "Others")],
    aes(
        x = windowCharge,
        y = IUPRED2,
        size = ifelse(
            JMJD6_substrate_flag == "JMJD6_substrate",
            oxK_ratio + 0.01,
            0.01
        ),
        fill = JMJD6_substrate_flag,
        alpha = JMJD6_substrate_flag,
        shape = JMJD6_substrate_flag
    )
) +
    geom_jitter(width = 0.04, height = 0) +
    theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "vertical"
    ) +
    scale_size_continuous(name = "Hydroxylation rate [%]") +
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

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

![](j3-properties-of-hydroxylation-sites_files/figure-html/analysis of sequence feature-1.png)<!-- -->

```r
ggplot(
    data = non.duplicated.stoichiometry.dt[order(JMJD6_substrate_flag != "Others")],
    aes(
        x = K_ratio_score,,
        y = IUPRED2,
        size = ifelse(
            JMJD6_substrate_flag == "JMJD6_substrate",
            oxK_ratio + 0.01,
            0.01
        ),
        fill = JMJD6_substrate_flag,
        alpha = JMJD6_substrate_flag,
        shape = JMJD6_substrate_flag
    )
) +
    geom_jitter(width = 0.04, height = 0) +
    theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "vertical"
    ) +
    scale_size_continuous(name = "Hydroxylation rate [%]") +
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
    xlab("K-score") +
    ylab("Disordedness") +
    ggtitle("Biophysical propety around K")
```

![](j3-properties-of-hydroxylation-sites_files/figure-html/analysis of sequence feature-2.png)<!-- -->



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
##  date     2021-12-08                  
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
