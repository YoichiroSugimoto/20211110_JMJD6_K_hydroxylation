---
title: "Significance of hydroxylation sites"
author: "Yoichiro Sugimoto"
date: "30 November, 2021"
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


# Load packages



```r
temp <- sapply(list.files("../functions", full.names = TRUE), source)
set.seed(0)

library("readxl")
library("exact2x2")
```

```
## Error in get(genname, envir = envir) : object 'testthat_print' not found
```

```r
j2.4.input.dir <- file.path("../../data/j2-PTM-stoichiometry/PTM-significance")

results.dir <- file.path("../../results")
j2.res.dir <- file.path(results.dir, "j2-PTM-stoichiometry")
j2.4.res.dir <- file.path(j2.res.dir, "j2-4-significance-of-PTM")

create.dirs(c(
    j2.4.res.dir
))
```

# Analysis



```r
count.dt <- read_xlsx(
    file.path(
        j2.4.input.dir,
        "Peptide_PD_data_for_fishers_test.xlsx"
    ),
    sheet = "new data"
) %>% as.data.table

setnames(
    count.dt,
    old = c(
        "(#K) Peptide PD: HeLa_wt", "(#K) Peptide PD: HeLa_J6KO",
        "(#Kox) Peptide PD: HeLa_wt", "(#Kox)Peptide PD: HeLa_J6KO"
    ),
    new = c("K_all_WT", "K_all_J6KO", "K_ox_WT", "K_ox_J6KO")
)

count.dt[, `:=`(
    K_unmodified_WT = K_all_WT - K_ox_WT,
    K_unmodified_J6KO = K_all_J6KO - K_ox_J6KO
)]

count.dt[
   ,
    c("odds_ratio", "fisher_p") := matrix(
        c(K_ox_WT, K_unmodified_WT, K_ox_J6KO, K_unmodified_J6KO),
        nrow = 2,
        dimnames = list(
            K_modification = c("hydroxylated", "unmodified"),
            JMJD6 = c("WT", "J6KO")
        )
    ) %>%
        {list(
             odds_ratio = fisher.test(.)$estimate,
             p_value = fisher.test(.)$p.value
         )},
    by = seq_len(nrow(count.dt))
]

count.dt[
   ,
    boschloo_p := boschloo(
        x1 = K_ox_WT, n1 = K_all_WT,
        x2 = K_ox_J6KO, n2 = K_all_J6KO,
        alternative = "two.sided"
    )$p.value,
    by = seq_len(nrow(count.dt))
]


count.dt[, `:=`(
    fisher_fdr = p.adjust(fisher_p, method = "fdr"),
    boschloo_fdr = p.adjust(boschloo_p, method = "fdr")
)]

fwrite(
    count.dt,
    file.path(j2.4.res.dir, "significance-of-hydroxylation-site-v2.csv")
)
```

