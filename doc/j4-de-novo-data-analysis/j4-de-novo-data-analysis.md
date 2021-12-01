---
title: "j4 de novo data analysis"
author: "Yoichiro Sugimoto"
date: "01 December, 2021"
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
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## 
## Attaching package: 'data.table'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
```

```r
source("../../R/j2-PTM-stoichiometry/functions/stoichiometry_functions.R")

results.dir <- file.path("../../results")
j4.res.dir <- file.path(results.dir, "j-d-de-novo")

create.dirs(c(
    j4.res.dir
))

processors <- 8
```

# Analysis


```r
j4.input.dir <- file.path("../../data/j4-de-novo-data")

peptide.wt.dt <- fread(file.path(j4.input.dir, "HeLa_peptide_de_novo_peptides.csv"))
peptide.j6ko.dt <- fread(file.path(j4.input.dir, "J6KO_peptide_de_novo_peptides.csv"))
jq1.wt.dt <- fread(file.path(j4.input.dir, "HeLa_JQ1_de_novo_peptides.csv"))
jq1.j6ko.dt <- fread(file.path(j4.input.dir, "J6KO_JQ1_de_novo_peptides.csv"))

peptide.wt.dt[, `:=`(PD = "peptide", JMJD6 = "WT")]
peptide.j6ko.dt[, `:=`(PD = "peptide", JMJD6 = "JMJD6KO")]
jq1.wt.dt[, `:=`(PD = "JQ1", JMJD6 = "WT")]
jq1.j6ko.dt[, `:=`(PD = "JQ1", JMJD6 = "JMJD6KO")]

all.dt <- rbindlist(list(peptide.wt.dt, peptide.j6ko.dt, jq1.wt.dt, jq1.j6ko.dt))
all.dt[
  , conf_score := list(str_split(get("local confidence (%)"), " ")), by = seq_len(nrow(all.dt))
]

all.dt[, `:=`(
    cell = paste0(PD, "_", JMJD6),
    Unique = "Y",
    PTM = "Oxidation (K)", # keep all KOH
    Accession = 1
)]

## Dummy data
all.protein.feature.per.pos.dt <- data.table(Accession = 1)

all.dt <- filterPeptideData(all.dt)
```

```
## Loading required package: Biostrings
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:data.table':
## 
##     first, second
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:data.table':
## 
##     shift
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## [1] "Input peptides per cells"
## cell
##     JQ1_JMJD6KO          JQ1_WT peptide_JMJD6KO      peptide_WT 
##          136947          128674          608508          585643 
## [1] "Peptides per cells after filteration 1"
## cell
##     JQ1_JMJD6KO          JQ1_WT peptide_JMJD6KO      peptide_WT 
##          136947          128674          608508          585643 
## [1] "Peptides per cells after filtration 2"
## cell
##     JQ1_JMJD6KO          JQ1_WT peptide_JMJD6KO      peptide_WT 
##           99960           96874          522394          500825 
## [1] "After the filtration of acetylation (filteration 3)"
## [1] "All K related modifications"
## 
##         K K(+15.99) K(+56.03) K(+72.02) 
##    431649     87034    917510     99962 
## [1] "All R related modifications"
## 
##       R 
## 1086699 
## [1] "All PTMs"
## 
## Oxidation (K) 
##       1220053 
## [1] "Peptides per cells after filtration 4"
## cell
##     JQ1_JMJD6KO          JQ1_WT peptide_JMJD6KO      peptide_WT 
##           98648           95262          513079          491497 
## [1] "Peptides per cells after filtration 5"
## cell
##     JQ1_JMJD6KO          JQ1_WT peptide_JMJD6KO      peptide_WT 
##           91144           89575          492429          471907 
## [1] "Filtration 5 (PTM assignment)"
## [1] "All K related modifications"
## 
##         K K(+15.99) K(+56.03) 
##    442686     59606    844950 
## [1] "All R related modifications"
## 
##       R 
## 1000409 
## [1] "All PTMs"
## 
## Oxidation (K) 
##       1145055
```

```r
## Peptides which does not contain K were excluded
all.dt <- all.dt[str_count(Peptide, "K") > 0]

all.dt[, `:=`(
  PSM_hydroxy = gsub("(K\\(\\+15.99\\)|K\\(\\+72.02\\))", "x", Peptide) %>%
      {gsub("\\s*\\([^\\)]+\\)","", .)} %>%
      {gsub("x","K(+15.99)", .)}
)]

fwrite(
    all.dt[
        grepl("K\\(\\+15.99\\)", PSM_hydroxy)
    ][JMJD6 == "WT"][
      , .(JMJD6, PSM_hydroxy)
    ],
    file.path(j4.res.dir, "de-novo-hydroxylation-site.txt"),
    sep = "\t"
)

plotOxDistribution <- function(all.dt, window.size, denovo.score.cutoff){
    
    all.dt <- all.dt[
        length >= window.size &
        get("Denovo Score") >= denovo.score.cutoff
    ]

    K.probs.list <- mclapply(
        AAStringSet(all.dt[, raw_peptide]),
        letterFrequencyInSlidingView,
        view.width = window.size, letters = "K", as.prob = TRUE,
        mc.cores = processors
    )

    k.max.probs <- mclapply(
        K.probs.list,
        max,
        mc.cores = processors
    ) %>% unlist

    centre.longest.stretch <- function(z){
        consec.dt <- with(rle(z), {
            ok <- values == TRUE
            ends <- cumsum(lengths)[ok]
            starts <- ends - lengths[ok] + 1
            cbind(starts, ends)
        }) %>% data.table

        ##consec.dt[ends - starts == max(consec.dt[, ends - starts])][1, round((starts + ends)/2)] ## This one return center
        consec.dt[ends - starts == max(consec.dt[, ends - starts])][1, ends] ## This one return the last occurence
    }

    k.max.first.index <- mcmapply(
        function(x, y){
            centre.longest.stretch(x[, "K"] == y)
        },
        x = K.probs.list,
        y = k.max.probs,
        mc.cores = processors
    )

    oxK.probs.list <- mclapply(
        AAStringSet(all.dt[, k_info_peptide]),
        letterFrequencyInSlidingView,
        view.width = window.size, letters = "H", as.prob = TRUE,
        mc.cores = processors
    )

    oxK.pos.list <- mclapply(
        AAStringSet(all.dt[, k_info_peptide]),
        letterFrequencyInSlidingView,
        view.width = 1, letters = "H", as.prob = TRUE,
        mc.cores = processors
    )

    ra.oxK.score <- mcmapply(
        function(x, y){
            ifelse(x[, "H"] == 1, as.integer(y), 0) %>%
                {zoo::rollapply(., window.size, sum)}
        },
        oxK.pos.list,
        all.dt[, conf_score],
        SIMPLIFY = FALSE,
        mc.cores = processors
    )

    ra.oxK.max.score <- mclapply(
        ra.oxK.score,
        max,
        mc.cores = processors
    ) %>% unlist

    ra.oxk.max.score.index <- mcmapply(
        function(x, y){
            centre.longest.stretch(x == y) 
        },
        ra.oxK.score,
        ra.oxK.max.score,
        mc.cores = processors
    )

    K.pos.list <- mclapply(
        AAStringSet(all.dt[, raw_peptide]),
        letterFrequencyInSlidingView,
        view.width = 1, letters = "K", as.prob = TRUE,
        mc.cores = processors
    )

    score.dt <- mcmapply(
        function(conf.score, oxk.score, oxk.index, kmax.index, k.pos, oxk.pos){
            shift.size <- window.size - 1
            
            if(oxk.score != 0){
                out <- case_when(
                    oxk.pos[, "H"][oxk.index:(oxk.index + shift.size)] == 1 ~
                        as.numeric(conf.score[oxk.index:(oxk.index + shift.size)]),
                    k.pos[, "K"][oxk.index:(oxk.index + shift.size)] == 1 ~ 1,
                    TRUE ~ 0
                )
            } else {
                out <- case_when(
                    k.pos[, "K"][kmax.index:(kmax.index + shift.size)] == 1 ~ 1,
                    TRUE ~ 0
                )
            }
            data.table(matrix(out, nrow = 1))
        },
        conf.score = all.dt[, conf_score],
        oxk.score = ra.oxK.max.score,
        oxk.index = ra.oxk.max.score.index,
        kmax.index = k.max.first.index,
        k.pos = K.pos.list,
        oxk.pos = oxK.pos.list,
        SIMPLIFY = FALSE,
        mc.cores = processors
    ) %>%
        rbindlist

    all.dt[, `:=`(
        max_K_ratio = k.max.probs,
        max_K_index = k.max.first.index,
        max_oxK_score = ra.oxK.max.score,
        max_oxK_index = ra.oxk.max.score.index
    )]

    all.score.dt <- cbind(all.dt, score.dt)

    all.score.file <- file.path("../../results", paste0(
                                                  "all_score_",
                                                  window.size, "_", denovo.score.cutoff
                                                , ".csv"
                                              ))
    
    fwrite(
        all.score.dt,
        all.score.file
    )

    all.score.dt <- fread(file.path(all.score.file))

    all.score.dt[, max_K_ratio2 := case_when(
                       max_K_ratio >= 0.8 ~ "[0.8, 1]",
                       max_K_ratio >= 0.6 ~ "[0.6, 0.8)",
                       max_K_ratio >= 0.4 ~ "[0.4, 0.6)",
                       max_K_ratio >= 0.2 ~ "[0.2, 0.4)",
                       max_K_ratio >= 0 ~ "[0, 0.2)"
                   )]

    print(all.score.dt[, .N, by = list(JMJD6, max_K_ratio2)])

    all.score.dt <- all.score.dt[
        order(- max_K_ratio2, - max_oxK_score)
    ][max_K_ratio > 0.1]

    all.score.dt[, sample_id := 1:.N, by = list(JMJD6)]

    m.all.score.dt <- melt(
        all.score.dt,
        id.vars = c("sample_id", "Peptide", "max_oxK_score", "max_K_ratio2", "JMJD6"),
        measure.vars = paste0("V", 1:window.size),
        value.name = "conf_score",
        variable.name = "position_in_n_mer"
    )

    m.all.score.dt[, `:=`(
        conf_score = ifelse(conf_score == 0, NA, conf_score),
        JMJD6 = factor(JMJD6, levels = c("WT", "JMJD6KO")),
        max_K_ratio2 = factor(max_K_ratio2, levels = c("[0.8, 1]", "[0.6, 0.8)", "[0.4, 0.6)", "[0.2, 0.4)", "[0, 0.2)"))
    )]

    m.all.score.dt[, `:=`(
        total_N_per_group = .N,
        rank_per_group = 1:.N
    ), by = list(JMJD6, max_K_ratio2, position_in_n_mer)]

    m.all.score.dt[, ratio_in_group := rank_per_group/total_N_per_group]

    plot.with.th <- function(ratio.th){

        g1 <- ggplot(
            m.all.score.dt[JMJD6 == "WT"][ratio_in_group < ratio.th],
            aes(
                x = position_in_n_mer,
                y = - sample_id,
                fill = conf_score
            )
        ) +
            geom_tile() +
            facet_grid(max_K_ratio2 ~ ., scales = "free") +
            scale_fill_gradient2(midpoint = 1, low = "white", mid = "white", high = "#000080", na.value = "white") +
            theme(
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.title = element_text()
            ) +
            xlab("Position in peptides") +
            ggtitle("WT")

        g2 <- ggplot(
            m.all.score.dt[JMJD6 == "JMJD6KO"][ratio_in_group < ratio.th],
            aes(
                x = position_in_n_mer,
                y = - sample_id,
                fill = conf_score
            )
        ) +
            geom_tile() +
            facet_grid(max_K_ratio2 ~ ., scales = "free") +
            scale_fill_gradient2(midpoint = 1, low = "white", mid = "white", high = "#000080", na.value = "white") +
            theme(
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.title = element_text()
            ) +
            xlab("Position in peptides") +
            ggtitle("JMJD6KO")

        library("cowplot")

        prow <- plot_grid(
            g1 + theme(legend.position="none"),
            g2 + theme(legend.position="none"),
            align = 'vh',
            hjust = -1,
            nrow = 1
        )

        legend <- get_legend(
            ## create some space to the left of the legend
            g1 + theme(legend.box.margin = margin(0, 0, 0, 12))
        )

        pa <- plot_grid(prow, legend, rel_widths = c(3, .4))

        title <- ggdraw() + 
            draw_label(
                paste0(
                    "Oxidized K distribution\n widow size: ",
                    window.size, "; De novo score >=: ", denovo.score.cutoff,
                    "; top ", round(ratio.th * 100), "% peptides shown"
                ),
                fontface = 'bold',
                x = 0,
                hjust = 0
            ) +
            theme(
                plot.margin = margin(0, 0, 0, 7)
            )

        plot_grid(
            title, pa,
            ncol = 1,
            rel_heights = c(0.1, 1)
        )

    }

    pg1 <- plot.with.th(1.2)
    pg2 <- plot.with.th(0.25)

    print(pg1)
    print(pg2)

    return()
}
```



```r
temp <- plotOxDistribution(
    all.dt = all.dt,
    window.size = 6,
    denovo.score.cutoff = 50
)
```

```
##       JMJD6 max_K_ratio2      N
##  1:      WT     [0.8, 1]    243
##  2:      WT   [0.6, 0.8)   3533
##  3:      WT   [0.4, 0.6)  20671
##  4:      WT   [0.2, 0.4)  98033
##  5:      WT     [0, 0.2) 261928
##  6: JMJD6KO     [0.8, 1]    281
##  7: JMJD6KO   [0.4, 0.6)  21229
##  8: JMJD6KO   [0.6, 0.8)   3505
##  9: JMJD6KO   [0.2, 0.4) 103023
## 10: JMJD6KO     [0, 0.2) 277327
```

```
## 
## ********************************************************
```

```
## Note: As of version 1.0.0, cowplot does not change the
```

```
##   default ggplot2 theme anymore. To recover the previous
```

```
##   behavior, execute:
##   theme_set(theme_cowplot())
```

```
## ********************************************************
```

![](j4-de-novo-data-analysis_files/figure-html/window size 6-1.png)<!-- -->![](j4-de-novo-data-analysis_files/figure-html/window size 6-2.png)<!-- -->



```r
temp <- plotOxDistribution(
    all.dt = all.dt,
    window.size = 10,
    denovo.score.cutoff = 50
)
```

```
##      JMJD6 max_K_ratio2      N
## 1:      WT   [0.4, 0.6)   9504
## 2:      WT   [0.2, 0.4) 114718
## 3:      WT     [0, 0.2) 168828
## 4:      WT   [0.6, 0.8)    104
## 5: JMJD6KO   [0.4, 0.6)   9098
## 6: JMJD6KO   [0.2, 0.4) 122315
## 7: JMJD6KO     [0, 0.2) 182058
## 8: JMJD6KO   [0.6, 0.8)    115
```

![](j4-de-novo-data-analysis_files/figure-html/windows size 10-1.png)<!-- -->![](j4-de-novo-data-analysis_files/figure-html/windows size 10-2.png)<!-- -->
