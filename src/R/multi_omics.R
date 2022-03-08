#===================================================================================================

#Multi-Omics time-course plot (Figure 6a)

#===================================================================================================
UI_DEMET_INTERVAL <- c("de00_48", "de00_54", "de00_60", "de00_66", "de00_72")
TIME_POINT <- c("de0h", "de48h", "de54h", "de60h", "de66h", "de72h")

#-------------------------------------------------
#GATA6 expression
#-------------------------------------------------
#loding the CAGE expression data
expMatrix <- read.table(
    "iPS_DE_diff_F5.TPM.txt",
    header = TRUE,
    row.name = NULL,
    sep = "\t",
    stringsAsFactors = FALSE
)

expMatrix <- expMatrix[, c(1, 3, 10, 17, 5, 12, 19, 8, 15, 22, 7, 14, 21, 2, 9, 16, 6, 13, 20)] # reorder based on time point

#gene
promoterID <- expMatrix[, "prmtrID"]
gene_name <- sapply(strsplit(promoterID, "@"), function(x){ x[2] })

#sum tpm by gene
expMatrix <- cbind(expMatrix, gene = gene_name)

expMatrix %>%
    group_by(gene) %>%
    select(-prmtrID) %>%
    summarise_each(dplyr::funs(sum)) -> sum_expMatrix

#Compute mean and sd of replicates
mean.sum_expMatrix <- cbind(
    apply(sum_expMatrix[, c(2:4)], 1, mean),
    apply(sum_expMatrix[, c(5:7)], 1, mean),
    apply(sum_expMatrix[, c(8:10)], 1, mean),
    apply(sum_expMatrix[, c(11:13)], 1, mean),
    apply(sum_expMatrix[, c(14:16)], 1, mean),
    apply(sum_expMatrix[, c(17:19)], 1, mean)
)
colnames(mean.sum_expMatrix) <- c(0, 48, 54, 60, 66, 72)
rownames(mean.sum_expMatrix) <- sum_expMatrix$gene

##compute sd of replicates
sd.sum_expMatrix <- cbind(
    apply(sum_expMatrix[, c(2:4)], 1, sd),
    apply(sum_expMatrix[, c(5:7)], 1, sd),
    apply(sum_expMatrix[, c(8:10)], 1, sd),
    apply(sum_expMatrix[, c(11:13)], 1, sd),
    apply(sum_expMatrix[, c(14:16)], 1, sd),
    apply(sum_expMatrix[, c(17:19)], 1, sd)
)
colnames(sd.sum_expMatrix) <- c(0, 48, 54, 60, 66, 72)
rownames(sd.sum_expMatrix) <- sum_expMatrix$gene

goi <- "GATA6"
goi.mean.sum_expMatrix <- mean.sum_expMatrix[rownames(mean.sum_expMatrix) %in% goi,, drop = FALSE]
goi.sd.sum_expMatrix <- sd.sum_expMatrix[rownames(sd.sum_expMatrix) %in% goi,, drop = FALSE]
#-------------------------------------------------
#M-values at demethylated regions
#-------------------------------------------------
extension <- c(0, 1)

controls <- c(1, 1, 1, 1, 1)
treatments <- c(2, 3, 4, 5, 6)

dmr.gr_list <- NULL
for(i in seq(length(controls))){
    bed <- dmr2bed(
        data_matrix = selDataMatrix,
        ControlColnum = controls[i],
        TreatmentColnum = treatments[i],
        MethylDemethyl = "Demethyl",
        extension = extension
    )
    gr <- bed2gr(bed)
    dmr.gr_list <- c(dmr.gr_list, list(gr))
}

names(dmr.gr_list) <- c("Demet.0h_48h", "Demet.0h_54h", "Demet.0h_60h", "Demet.0h_66h", "Demet.0h_72h")

dmr.probeID <- lapply(dmr.gr_list, function(x){ x$id })

#uninherited
uninherit_demet_grl <- GRangesList(
    Demet.0h_48h = dmr.gr_list[[1]],
    Demet.0h_54h = dmr.gr_list[[2]][!(dmr.probeID[[2]] %in% unlist(dmr.probeID[1]))],
    Demet.0h_60h = dmr.gr_list[[3]][!(dmr.probeID[[3]] %in% unlist(dmr.probeID[1:2]))],
    Demet.0h_66h = dmr.gr_list[[4]][!(dmr.probeID[[4]] %in% unlist(dmr.probeID[1:3]))],
    Demet.0h_72h = dmr.gr_list[[5]][!(dmr.probeID[[5]] %in% unlist(dmr.probeID[1:4]))]
)

#-------------------------------------------------
#transcript at demethylated regions
#-------------------------------------------------
#MDR expression table list
data_dir <- "dataMethyl850_DE/exp_table/"
files <- c(
    "iPS-DE_0.48_demethyl.TPM.txt",
    "iPS-DE_0.54_demethyl.TPM.txt",
    "iPS-DE_0.60_demethyl.TPM.txt",
    "iPS-DE_0.66_demethyl.TPM.txt",
    "iPS-DE_0.72_demethyl.TPM.txt"
)

group_ids <- c("iPS_DE_00h","iPS_DE_48h","iPS_DE_54h","iPS_DE_60h","iPS_DE_66h","iPS_DE_72h")

#load CAGE expression tables
exp_list <- NULL
for(i in seq(length(files))){
    #load a expression table
    exp_data <- read.table(
        file = paste0(data_dir, files[i]),
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    
    #Set header to sample name
    header2 <- NULL
    header <- colnames(exp_data)
    for(j in header){
        index <- SAMPLE_ID[, 1] %in% j
        if(!any(index)){
            sample_name <- "probeID"
        }else{
            sample_name <- SAMPLE_ID[SAMPLE_ID[, 1] %in% j, 2]
        }
        header2 <- c(header2, sample_name)
    }
    
    colnames(exp_data) <- header2
    
    #column reorder
    id_index <- 1
    reorder_index <- order(colnames(exp_data)[-id_index])
    
    exp_data2 <- cbind(probeID = exp_data[, id_index], exp_data[, - id_index][, reorder_index])
    
    #output list
    exp_list <- c(exp_list, list(exp_data2))
}
names(exp_list) <- UI_DEMET_INTERVAL

#mean and sd of rxpression
mean_expMatrix_list <- NULL
sd_expMatrix_list <- NULL
for(i in seq(length(exp_list))){
    exp_table <- exp_list[[i]]
    ##compute mean of replicates
    mean_expMatrix <- NULL
    sd_expMatrix <- NULL
    
    for(j in seq(length(group_ids))){
        #mean
        mean_exp <- select(exp_table, starts_with(group_ids[j])) %>% rowMeans()
        mean_expMatrix <- cbind(mean_expMatrix, mean_exp)
        #sd
        sd_exp <- select(exp_table, starts_with(group_ids[j])) %>% apply(., 1, sd)
        sd_expMatrix <- cbind(sd_expMatrix, sd_exp)
    }
    
    colnames(mean_expMatrix) <- paste0(TIME_POINT, "exp_mean")
    rownames(mean_expMatrix) <- exp_table$probeID
    mean_expMatrix_list <- c(mean_expMatrix_list, list(mean_expMatrix))
    
    colnames(sd_expMatrix) <- paste0(TIME_POINT, "exp_sd")
    rownames(sd_expMatrix) <- exp_table$probeID
    sd_expMatrix_list <- c(sd_expMatrix_list, list(sd_expMatrix))
}
names(mean_expMatrix_list) <- UI_DEMET_INTERVAL
names(sd_expMatrix_list) <- UI_DEMET_INTERVAL

#-------------------------------------------------#
#ATAC scores at demethylated regions
#-------------------------------------------------#
#load the bigwig data
atac.bws <- c("atac_merged_0h.bw", "atac_merged_48h.bw", "atac_merged_54h.bw", "atac_merged_60h.bw", "atac_merged_66h.bw", "atac_merged_72h.bw")

#compute the mean score at demethylated regions
atac_scores_list <- NULL
for(i in seq(length(uninherit_demet_grl))){
    demet_ext <- 250
    #Scores of atac
    atac_scores <- data.frame(id = uninherit_demet_grl[[i]]$id, mval = score(uninherit_demet_grl[[i]]))
    for(j in seq(length(TIME_POINT))){
        #ATAC-seq
        atac.gr <- import.bw(
            con = atac.bws[j],
            format = "bigWig",
            as = "GRanges"
        )
        
        overlap_atac_index <- findOverlaps(uninherit_demet_grl[[i]] + demet_ext, atac.gr)
        overlap_atac.gr <- atac.gr[subjectHits(overlap_atac_index)]
        overlap_atac.gr$id <- uninherit_demet_grl[[i]][queryHits(overlap_atac_index)]$id
        
        id_score_atac <- data.frame(id = overlap_atac.gr$id, score = score(overlap_atac.gr))
        id_score_atac %>%
            group_by(id, add = FALSE) %>%
            summarise(mean = mean(score)) %>%
            left_join(atac_scores, ., by = "id") -> atac_scores
        
        rm(atac.gr)
        rm(overlap_atac.gr)
        gc(gc())
    }
    names(atac_scores) <- c("id", "Mval", TIME_POINT)
    atac_scores_list <- c(atac_scores_list, list(atac_scores))
}
#-------------------------------------------------#
#ChIPmentation at demethylated regions
#-------------------------------------------------#
#load the bigwig data
chip.bws <- c("GATA6_0h.sorted.bw", "GATA6_48h.sorted.bw", "GATA6_54h.sorted.bw", "GATA6_60h.sorted.bw", "GATA6_66h.sorted.bw", "GATA6_72h.sorted.bw")


#compute the mean score at demethylated regions
chip_scores_list <- NULL
for(i in seq(length(uninherit_demet_grl))){
    demet_ext <- 250
    chip_scores <- data.frame(id = uninherit_demet_grl[[i]]$id, mval = score(uninherit_demet_grl[[i]]))
    for(j in seq(length(TIME_POINT))){
        #ChIPmentation
        chip.gr <- import.bw(
            con = chip.bws[j],
            format = "bigWig",
            as = "GRanges"
        )
        
        overlap_chip_index <- findOverlaps(uninherit_demet_grl[[i]] + demet_ext, chip.gr)
        overlap_chip.gr <- chip.gr[subjectHits(overlap_chip_index)]
        overlap_chip.gr$id <- uninherit_demet_grl[[i]][queryHits(overlap_chip_index)]$id
        
        id_score_chip <- data.frame(id = overlap_chip.gr$id, score = score(overlap_chip.gr))
        id_score_chip %>%
            group_by(id, add = FALSE) %>%
            summarise(mean = mean(score)) %>%
            left_join(chip_scores, ., by = "id") -> chip_scores
        
        rm(chip.gr)
        rm(overlap_chip.gr)
        gc(gc())
    }
    names(chip_scores) <- c("id", "Mval", TIME_POINT)
    chip_scores_list <- c(chip_scores_list, list(chip_scores))
}

#-------------------------------------------------#
#Integration of all omics data
#-------------------------------------------------#
for(i in seq(length(uninherit_demet_grl))){
    cat(paste0("\t",names(uninherit_demet_grl)[i], "\n"))
    # M value processing
    met.mval <- selDataMatrix[rownames(selDataMatrix) %in% uninherit_demet_grl[[i]]$id, 1:6]
    
    nega_delta.mval <- data.frame(
        met.mval[1] - met.mval[1],
        met.mval[1] - met.mval[2],
        met.mval[1] - met.mval[3],
        met.mval[1] - met.mval[4],
        met.mval[1] - met.mval[5],
        met.mval[1] - met.mval[6]
    )
    
    names(nega_delta.mval) <- TIME_POINT
    
    mval_mean <- sapply(nega_delta.mval, mean)
    
    # CAGE processing
    cage.tpm <- as.data.frame(mean_expMatrix_list[[i]])
    cage.tpm[cage.tpm == 0] <- 0.001
    
    cage_lfc <- data.frame(
        de0h = log(cage.tpm[1] / cage.tpm[1], 2),
        de48h = log(cage.tpm[2] / cage.tpm[1], 2),
        de54h = log(cage.tpm[3] / cage.tpm[1], 2),
        de60h = log(cage.tpm[4] / cage.tpm[1], 2),
        de66h = log(cage.tpm[5] / cage.tpm[1], 2),
        de72h = log(cage.tpm[6] / cage.tpm[1], 2)
    )
    names(cage_lfc) <- TIME_POINT
    
    cage_mean <- sapply(cage_lfc, mean)
    
    # atac processing
    atac_scores2 <- tidyr::replace_na(atac_scores_list[[i]], replace = list(de0h = 0.001, de48h = 0.001, de54h = 0.001, de60h = 0.001, de66h = 0.001, de72h = 0.001))
    
    atac_lfc <- data.frame(
        id = atac_scores2$id,
        de0h = log(atac_scores2$de0h / atac_scores2$de0h, 2),
        de48h = log(atac_scores2$de48h / atac_scores2$de0h, 2),
        de54h = log(atac_scores2$de54h / atac_scores2$de0h, 2),
        de60h = log(atac_scores2$de60h / atac_scores2$de0h, 2),
        de66h = log(atac_scores2$de66h / atac_scores2$de0h, 2),
        de72h = log(atac_scores2$de72h / atac_scores2$de0h, 2)
    )
    
    atac_mean <- sapply(atac_lfc[-1], mean)
    
    # chip processing
    chip_scores2 <- tidyr::replace_na(chip_scores_list[[i]], replace = list(de0h = 0.001, de48h = 0.001, de54h = 0.001, de60h = 0.001, de66h = 0.001, de72h = 0.001))
    
    chip_lfc <- data.frame(
        id = chip_scores2$id,
        de0h = log(chip_scores2$de0h / chip_scores2$de0, 2),
        de48h = log(chip_scores2$de48h / chip_scores2$de0h, 2),
        de54h = log(chip_scores2$de54h / chip_scores2$de0h, 2),
        de60h = log(chip_scores2$de60h / chip_scores2$de0h, 2),
        de66h = log(chip_scores2$de66h / chip_scores2$de0h, 2),
        de72h = log(chip_scores2$de72h / chip_scores2$de0h, 2)
    )
    
    chip_mean <- sapply(chip_lfc[-1], mean)
    
    # Line plot
    library(tidyverse)
    #data integration
    second_rate = 0.01 
    df <- as.data.frame(rbind(goi.mean.sum_expMatrix * second_rate, atac_mean, chip_mean, cage_mean, mval_mean))
        names(df) <- c(0, 48, 54, 60, 66, 72)
    df$data <- factor(c("GATA6", "ATAC", "ChIP", "Transcript", "M-value"), levels = c("GATA6", "ATAC", "ChIP", "Transcript", "M-value"))
    
    df2 <- pivot_longer(df, cols = -data, values_to = "value", names_to = "hours")
    df2$hours <- as.numeric(df2$hours)

    library(ggplot2)
    library(ggsci)
    g <- ggplot(df2, aes(x = hours, y = value, color = data, group = data))
    g <- g + geom_line(size = 1.5)
    g <- g + geom_point(size = 2)
    g <- g + scale_y_continuous(sec.axis = sec_axis(~. / second_rate, name = "GATA6 (TPM)"))
    g <- g + labs(x = "Time Point(h)", y = "ATAC, ChIP, TR, M-value(og2FC, vs. 0 h)")
    g <- g +    theme(
            axis.text = element_text(color = "black", size = 20),    # axis texit setting
            axis.ticks = element_line(color = "black", size = .5),    # Setting of scale ticks line
            axis.title = element_text(color = "black", size = 20),    # Setting of scale ticks texit
            panel.background = element_rect(fill = "white", color = NA),    # panel background color
            panel.border = element_rect(fill = NA, color = "black"),    # panel border
            panel.grid.major = element_line(color = "gray90"),    # panel major grid
            panel.grid.minor = element_line(color = "gray90"), # panel minor grid
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_rect(fill = NA),
            legend.text = element_text(size=20)
        )
    g <- g + scale_color_jama(alpha=c(0.3, 1, 1, 1, 1))
    g
    
    ggsave(
        g,
        file = paste0(names(uninherit_demet_grl)[i], "omics_plot.pdf"),
        width = 9,
        heigh = 7
    )
}