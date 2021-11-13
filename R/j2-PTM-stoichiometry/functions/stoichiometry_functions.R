countPTM <- function(all.pp.dt, peptide.col = "Peptide"){

    all.dt <- copy(all.pp.dt)
    
    print("All K related modifications")
    all.letters <- all.dt[, unlist(str_split(get(peptide.col), "(?=[[:alpha:]])"))]
    all.letters <- all.letters[all.letters != ""]
    print(table(all.letters[grepl("^K", all.letters)]))
    all.dt[, aa_list := str_split(Peptide, "(?=[[:alpha:]])")]

    print("All R related modifications")
    all.letters <- all.dt[, unlist(str_split(get(peptide.col), "(?=[[:alpha:]])"))]
    all.letters <- all.letters[all.letters != ""]
    print(table(all.letters[grepl("^R", all.letters)]))
    all.dt[, aa_list := str_split(Peptide, "(?=[[:alpha:]])")]

    print("All PTMs")
    PTMs <- all.dt[, unlist(str_split(PTM, "; "))]
    PTMs <- PTMs[PTMs != ""]
    print(table(sort(PTMs)))

    return()
}


filterPeptideData <- function(all.pp.dt){
    require("Biostrings")

    all.dt <- copy(all.pp.dt)

    print("Input peptides per cells")
    print(all.dt[, table(cell)])

    ## F1.
    all.dt <- all.dt[Unique == "Y"]
    print("Peptides per cells after filteration 1")
    print(all.dt[, table(cell)])    

    ## F2
    all.dt <- all.dt[
        Area > 0 &
        Accession %in% all.protein.feature.per.pos.dt[, unique(Accession)]
    ]

    print("Peptides per cells after filtration 2")
    print(all.dt[, table(cell)])

    ## F3. Remove Acetylation from modification
    all.dt[, corrected_peptide := gsub("\\(\\+42\\.01\\)", "", Peptide)]
    print("After the filtration of acetylation (filteration 3)")
    temp <- countPTM(all.dt, peptide.col = "corrected_peptide")

    ## F4. Filter out peptides with C-term propionylation and peptides which do not have C-term K, KOH or R
    all.dt <- all.dt[
    (
        str_count(corrected_peptide, "K\\(\\+56\\.03\\)$") +
        str_count(corrected_peptide, "K\\(\\+72\\.02\\)$") +
        str_count(corrected_peptide, "R\\(\\+28\\.03\\)$") == 0
    ) &
    (
        str_count(corrected_peptide, "(R$|R\\(\\+14\\.02\\)$)") == 1 |
        str_count(corrected_peptide, "(K$|K\\(\\+15\\.99\\)$)") == 1 
    )
    ]

    print("Peptides per cells after filtration 4")
    print(all.dt[, table(cell)])

    ## F5.
    ##  Allow only 1 misclevage by trypsin 
    ## - At most one occurrence of K, K(+15.99), R, R(+14.02) in the peptide
    ## - Note that the following are allowed in the peptides
    ##   - C terminal K, R, KP, RP
    ##   - Propionylated K: K(+56.03), K(+72.02)
    ##   - Dimethyl R: R(+28.03)

    all.dt[
      , Peptide_for_filtration := gsub(
            "(K$|K\\(\\+15\\.99\\)$|R$|R\\(\\+14\\.02\\)$)",
            "X",
            corrected_peptide
            ## C-term K or R is OK
        ) %>%        
            {gsub("K\\(\\+56\\.03\\)", "X", .)} %>%
            {gsub("K\\(\\+72\\.02\\)", "X", .)} %>%
            {gsub("R\\(\\+28\\.03\\)$", "X", .)} %>%
            ## The above modification prevents clevage by trypsin
            {gsub("\\s*\\([^\\)]+\\)", "", .)} %>%
            {gsub("KP", "XZ", .)} %>%
            {gsub("RP", "XZ", .)}
    ]

    all.dt <- all.dt[str_count(Peptide_for_filtration, "(K|R)") < 2]
    all.dt[, Peptide_for_filtration := NULL]

    print("Peptides per cells after filtration 5")
    print(all.dt[, table(cell)])

    all.dt[, PTMs := str_split(PTM, "; ")]

    all.dt <- all.dt[
      , corrected_peptide := case_when(
            !any(PTMs[[1]] == "Oxidation (K)") ~
                gsub("K\\(\\+15.99\\)","K", corrected_peptide),
            TRUE ~ corrected_peptide
        ) %>% {
            case_when(
                !any(PTMs[[1]] == "Oxidised Propionylation") ~
                    gsub("K\\(\\+72.02\\)","K", .),
                TRUE ~ .
            )} %>% {
                case_when(
                    !any(PTMs[[1]] == "Methylation(R)") ~
                        gsub("R\\(\\+14\\.02\\)", "R", .),
                    TRUE ~ .
                )
            } %>% {
                case_when(
                    !any(PTMs[[1]] == "Dimethylation(R)") ~
                        gsub("R\\(\\+28\\.03\\)", "R", .),
                    TRUE ~ .
                )},
        by = seq_len(nrow(all.dt))
    ]

    print("Filtration 5 (PTM assignment)")
    temp <- countPTM(all.dt, peptide.col = "corrected_peptide")

    ## For calculation of stoichiometry
    all.dt[, `:=`(
        raw_peptide = gsub("\\s*\\([^\\)]+\\)","", Peptide),
        k_info_peptide =
            gsub(
                names(AMINO_ACID_CODE) %>%
                {.[. != "K"]} %>%
                {paste(., collapse = "|")} %>%
                {paste0("(", ., ")")}, 
                "A", corrected_peptide # Non K amino acids will be replaced with A
            ) %>%
            {gsub("K\\(\\+56.03\\)", "K", .)} %>% 
            {gsub("(K\\(\\+15.99\\)|K\\(\\+72.02\\))", "H", .)} %>% # H: hydroxylared K
            {gsub("\\s*\\([^\\)]+\\)","", .)},
        r_info_peptide =
            gsub(
                names(AMINO_ACID_CODE) %>%
                {.[. != "R"]} %>%
                {paste(., collapse = "|")} %>%
                {paste0("(", ., ")")}, 
                "A", corrected_peptide # Non R amino acids will be replaced with A
            ) %>%
            {gsub("R\\(\\+14\\.02\\)", "M", .)} %>% # M: Monoaceyltation 
            {gsub("R\\(\\+28\\.03\\)", "D", .)} %>% # D: dimethylation
            {gsub("\\s*\\([^\\)]+\\)","", .)}
    )]

    all.dt[, `:=`(
        K_index = str_locate_all(k_info_peptide, "(K|H)") %>%
            {lapply(., function(x) x[,2])},
        oxK_index = str_locate_all(k_info_peptide, "H") %>%
            {lapply(., function(x) x[,2])},
        all_aa_index = str_locate_all(k_info_peptide, "(K|H|A)") %>%
            {lapply(., function(x) x[,2])},
        R_index = str_locate_all(r_info_peptide, "(R|M|D)") %>%
            {lapply(., function(x) x[,2])},
        mR_index = str_locate_all(r_info_peptide, "M") %>%
            {lapply(., function(x) x[,2])},
        dmR_index = str_locate_all(r_info_peptide, "D") %>%
            {lapply(., function(x) x[,2])}
        )]

    return(all.dt)
}

indexToPosition <- function(index.name, dt){
    
    dt <- copy(dt)
    dt[, peptide_id := 1:.N]
    dt[
      , n_index := length(get(paste0(index.name, "_index"))[[1]]),
        by = seq_len(nrow(dt)) 
    ]

    if(dt[, sum(n_index)] > 0){
        e.dt <- dt[rep(1:nrow(dt), times = n_index)]
        e.dt[, index_pos_id := 1:.N, by = peptide_id]
        e.dt[
          , position_in_peptide :=
                get(paste0(index.name, "_index"))[[1]][index_pos_id],
            by = seq_len(nrow(e.dt))
        ]
        e.dt[, `:=`(
            position = Start + position_in_peptide - 1,
            aa_type = index.name
        )]
        e.dt <- e.dt[, c(
            "cell", "Accession", "aa_type", "position", "Area", "n_feature"
        ), with = FALSE]

        pos.dt <- e.dt[, list(
            total_count = .N,
            total_area = sum(Area),
            total_n_feature = sum(n_feature)
        ), by = list(
               cell, 
               Accession, aa_type, position
           )]
    } else {
        pos.dt <- data.table()
    }

    return(pos.dt)
}
