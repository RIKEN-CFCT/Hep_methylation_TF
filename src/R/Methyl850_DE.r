#===================================================================================================

#preprocessing 

#===================================================================================================
library("InfiniumDiffMetMotR")

selDataMatrix <- lumiMethyNorm(
    fileName = "TableControl.txt",
    inputtype = "signal",
    sample_names = c("ID00", "ID48", "ID54", "D60", "ID66", "ID72", "ID96")
)

#===================================================================================================

# DMR b2b plot (Fig. 4B)

#===================================================================================================
#Count DMP
comp <- data.frame(
    DE0h_DE48h = c(1, 2),
    DE48h_DE54h = c(2, 3),
    DE54h_DE60h = c(3, 4),
    DE60h_DE66h = c(4, 5),
    DE66h_DE72h = c(5, 6)
)

M_ndiffprobes <- sapply(comp, function(x){
    sum(selDataMatrix[, as.numeric(x[1])] - selDataMatrix[, x[2]] <= -2)
})
D_ndiffprobes <- sapply(comp, function(x){
    sum(selDataMatrix[, as.numeric(x[1])] - selDataMatrix[, x[2]] >= 2)
})

nprobes_df <- data.frame(
    Time = c("0h-48h", "48h-54h", "54h-60h", "60h-66h", "66h-72h"),
    Methylated = M_ndiffprobes,
    Demethylated = D_ndiffprobes
)

#plot
nprobes_b2b_df <- tidyr::pivot_longer(nprobes_df, col = -Time,, names_to = "MetDemet", values_to = "nprobes")

library(dplyr)
the_order <- rev(as.character(unlist(
    nprobes_b2b_df %>%
        filter(MetDemet == "Methylated") %>%
        select(Time)
)))

theme <- theme(
    panel.background = element_blank(), # initialization
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "line"),
    panel.grid.major.x = element_line(linetype = "dotted", size = 0.3, color = "#3A3F4A")
)

g <- ggplot(nprobes_b2b_df, aes(x = Time))
g <- g + scale_x_discrete(limits = the_order)

#Methylated plot (left)
g.M <- g + geom_bar(
    data = nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]] == "Methylated",],
    aes(y = nprobes), fill = rgb(0.8, 0, 0.2, alpha = 0.9),
    stat = "identity"
)
g.M <- g.M + scale_y_reverse(limits = c(max(nprobes_b2b_df[["nprobes"]]) * 1.2, 0))
g.M <- g.M + ggtitle('Methylated')
g.M <- g.M + ylab("Number of probes")
g.M <- g.M + coord_flip()
g.M <- g.M + theme
g.M <- g.M + theme(
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_line(size = 0.5),
    axis.line.x = element_line(size = 0.5)
)

#Demethylated plot (Right)
g.D <- g + geom_bar(
    data = nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]] == "Demethylated",],
    aes(y = nprobes), fill = rgb(0, 0.7, 0, alpha = 0.9),
    stat = "identity"
)
g.D <- g.D + scale_y_continuous(limits = c(0, max(nprobes_b2b_df[["nprobes"]]) * 1.2))
g.D <- g.D + ggtitle('Demethylated')
g.D <- g.D + ylab("Number of probes")
g.D <- g.D + coord_flip()
g.D <- g.D + theme
g.D <- g.D + theme(
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_line(size = 0.5),
    axis.line.x = element_line(size = 0.5)
)

#labels (Center)
g.time <- g + geom_bar(
    data = nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]] == "Demethylated",],
    aes(y = 0, fill = 'white'),
    alpha = 0, stat = "identity"
)
g.time <- g.time + geom_text(aes(y = 0, label = as.character(Time)), size = 4)
g.time <- g.time + ggtitle("Time")
g.time <- g.time + coord_flip()
g.time <- g.time + ylab("")
g.time <- g.time + theme
g.time <- g.time + theme(panel.grid.major.x = element_blank())

g_all <- grid.arrange(
    g.M, g.time, g.D,
    ncol = 3,
    widths = c(4.25 / 10, 1.5 / 10, 4.25 / 10)
)

ggsave(g_all, file = "nDMP_b2bplot.pdf", width = 7, height = 7)

#===================================================================================================

#enrichment plot (Fig. S5B and S5C)

#===================================================================================================
#motif database construction
motif_names <- "GATA6"
motifDBList <- IMAGE_PWMlist
motif <- motifDBList[grep(motif_names, names(motifDBList))]

# parameter setting
outname <- "DE" #output file name
MethylDemethyl <- "Demethyl" #Indication of methylated regions or demethylated regions
COMP <- data.frame(DE_00_48 = c(1, 2), DE_48_54 = c(2, 3), DE_54_60 = c(3, 4), DE_60_66 = c(4, 5), , DE_66_72 = c(5, 6), DE_00_72 = c(1, 6))

#enrichment score computation
All_enrichments <- as.list(NULL) #List file to be input all nerichment scores
comp_names <- NULL #vector to be input comparison names

for(i in seq(ncol(COMP))){
    ControlColnum <- as.numeric(COMP[1, i])
    TreatmentColnum <- as.numeric(COMP[2, i])
    comp_names <- c(
        comp_names,
        paste0(colnames(selDataMatrix)[ControlColnum], "_vs_", colnames(selDataMatrix)[TreatmentColnum])
    )

    #extraction of differentially methylated probe IDs
    DMP_IDs <- DmpId(selDataMatrix = selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, MethylDemethyl = Demethyl)
    nDMP_IDs <- length(DMP_IDs)

    #extraction of DMP positions and stratified sampling
    probe_annotation <- EPICanno
    target_position <- na.omit(probeID2position(probe_IDs = DMP_IDs, anno_info = probe_annotation)) #conversion of DMP IDs to position
    randomProbe_IDs <- stratSampling(target_IDs = DMP_IDs, anno_info = probe_annotation) #stratified sampling for negative control
    random_position <- na.omit(probeID2position(probe_IDs = randomProbe_IDs, anno_info = probe_annotation)) #conversion of NC probe IDs to position
    positionsList <- list("target" = target_position, "random" = random_position) #integrate DMP positoins and NC positions

    #Sequence extraction
    library("BSgenome.Hsapiens.UCSC.hg19")
    tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep = ":"))
    genome <- eval(parse(text = tmp))

    seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
    sequences <- lapply(positionsList, function(x){ seqExtract(positions = x, genome = genome, seq_range) })

    #Write splitted sequences
    tempDir <- paste(outname, "_temp", sep = "")
    dir.create(tempDir)

    seqs <- sequences$target
    target_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "target")
    seqs <- sequences$random
    random_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "random")
    rm(sequences)
    rm(seqs)
    invisible(replicate(3, gc()))

    #motif search
    target_positionsList <- splitSeqMotDist(filenames = target_all_filenames, motif_list = motif, min.score = min.score)
    file.remove(target_all_filenames)
    gc()
    random_positionsList <- splitSeqMotDist(filenames = random_all_filenames, motif_list = motif, min.score = min.score)
    file.remove(random_all_filenames)
    gc()

    #computation of enrichment score
    target_mot_posi <- unlist(target_positionsList)
    ctrl_mot_posi <- unlist(random_positionsList)
    enrichment_scores <- enrichScoreDist(target_mot_posi, ctrl_mot_posi, seq_range = seq_range, motif_name = motif_names, nDMP_IDs = nDMP_IDs, outname = outname, plot_draw = FALSE)
    
    All_enrichments <- c(All_enrichments, list(enrichment_scores))

    file.remove(tempDir)
}

#Vidualization
ranks <- seq(seq_range[1], seq_range[2], length = 1001) # x-axis values
data <- sapply(All_enrichments, cbind)
rownames(data) <- ranks
colnames(data) <- comp_names

data.df <- as.data.frame(data)
data.df %>%
    mutate(Position = rownames(data.df)) %>%
    mutate_at(vars("Position"), as.integer) -> data.df

data.df2 <- tidyr::pivot_longer(data.df, cols = -Position, names_to = "comparison", values_to = "ES")

data_names <- c("ID00_vs_ID48", "ID48_vs_ID54", "ID54_vs_ID60", "ID60_vs_ID66", "ID66_vs_ID72", , "ID00_vs_ID72")
cols <- brewer.pal(7, "Dark2")

g_list <- NULL
for(i in seq(length(data_names))){
    data.df3 <- filter(data.df2, comparison == eval(data_names[i]))

    g <- ggplot(data.df3, aes(x = Position, y = ES, color = eval(cols[i])))
    g <- g + geom_line()
    g <- g + ylim(-0.1, 0.25)
    g <- g + xlab("Distance from CpG") # X-axis label
    g <- g + ylab("Enrichment Score") # Y-axis label
    g <- g + geom_hline(yintercept = 0) # wtire a horizontal line at intercept=0
    g <- g + geom_vline(xintercept = 0, color = "Gray51") # wtire a horizontal line at intercept=0
    g <- g + scale_colour_manual(values = c(cols[i]))
    ## Settings of axes
    g <- g + theme(
        axis.text.x = element_text(color = "black", size = 14), # X-axis texit setting
        axis.text.y = element_text(color = "black", size = 14), # Y-axis texit setting
        axis.ticks = element_line(color = "black", size = .5), # Setting of scale ticks line
        axis.title = element_text(color = "black", size = 14), # Setting of scale ticks texit
        panel.background = element_rect(fill = "white", color = NA), # panel background color
        panel.border = element_rect(fill = NA, color = "black"), # panel border
        panel.grid.major = element_blank(), # panel major grid
        panel.grid.minor = element_blank() # panel minor grid
        )
    ## Settings of legend
    g <- g + theme(legend.position = 'none')
    g_list <- c(g_list, list(g))
}

g_all <- plot_grid(plotlist = g_list, nrow = 1)
pdfout <- paste(outname, "_enrichment_plot.pdf", sep = "")
ggsave(g_all, file = pdfout, width = 42, height = 7)

#===================================================================================================

#overlap of demethylated probes(upset plot)(Fig. 4C)

#===================================================================================================
bed_list <- NULL
for( i in 2:6){
    bed <- dmr2bed(data_matrix = selDataMatrix, ControlColnum = 1, TreatmentColnum = i)
    bed_list <- c(bed_list, list(bed))
}
data <- lapply(bed_list, function(x){ x[, 4] })
names(data) <- c("Demet.0h_48h", "Demet.0h_54h", "Demet.0h_60h", "Demet.0h_66h", "Demet.0h_72h")

#upset plot
pdf("upset.pdf", width = 17.5, height = 7)
upset(fromList(data),
    sets = rev(names(data)),
    nsets = 5,
    order.by = "freq",
    set_size.show = TRUE,
    keep.order = TRUE,
    set_size.scale_max = 1000,
    sets.bar.color = rgb(0, 0.7, 0, alpha = 0.9)
)
dev.off()

#===================================================================================================
#
#methylation heatmap for uninheited demet probes (Fig. S5D)
#
#===================================================================================================
df <- selDataMatrix[rownames(selDataMatrix) %in% data[5], 1:6]
names(df) <- c("0 h", "48 h", "54 h", "60 h", "66 h", "72 h")

#heatmap
breaks <- seq(-2, 2, length = 500) 
mycols <- colorRamp2(breaks, colorRampPalette(brewer.pal(9, "PuBu"))(500))
ht <- tcHeatmap(df, mycols = mycols, z = TRUE)

pdf("DE0_72_M_heatmap.pdf", width = 7, height = 10.5)
draw(ht)
dev.off()
#===================================================================================================

#enrichment plot for uninherited demethylated regions (Fig. 4D)

#===================================================================================================
#uninherited
uninherit_demet_list <- list(
    Demet.0h_48h = data[[1]],
    Demet.0h_54h = data[[2]][!(data[[2]] %in% unlist(data[1]))],
    Demet.0h_60h = data[[3]][!(data[[3]] %in% unlist(data[1:2]))],
    Demet.0h_66h = data[[4]][!(data[[4]] %in% unlist(data[1:3]))],
    Demet.0h_72h = data[[5]][!(data[[5]] %in% unlist(data[1:4]))],
)


outname <- "uninherited" #output file name

#motif pwm
motif_names <- "GATA6"
motifDBList <- IMAGE_PWMlist
motif <- motifDBList[grep(motif_names, names(motifDBList))]

#enrichment score computation
All_enrichments <- as.list(NULL) #List file to be input all nerichment scores
comp_names <- NULL #vector to be input comparison names
for(i in seq(length(uninherit_demet_list))){
    # loop of the computation by comarison
    comp_names <- c(
        comp_names,
        gsub("ui", "", names(uninherit_demet_list)[i])
    )

    #extraction of differentially methylated probe IDs
    DMP_IDs <- uninherit_demet_list[[i]]
    nDMP_IDs <- length(DMP_IDs)
    allProbe_IDs <- rownames(selDataMatrix)

    #Extraction of DMP positions and stratified sampling
    probe_annotation <- InfiniumDiffMetMotR::EPICanno

    target_position <- na.omit(probeID2position(probe_IDs = DMP_IDs, anno_info = probe_annotation)) #conversion of DMP IDs to position
    randomProbe_IDs <- stratSampling(target_IDs = DMP_IDs, anno_info = probe_annotation) #stratified sampling for negative control
    random_position <- na.omit(probeID2position(probe_IDs = randomProbe_IDs, anno_info = probe_annotation)) #conversion of NC probe IDs to position
    positionsList <- list("target" = target_position, "random" = random_position) #integrate DMP positoins and NC positions

    #Sequence extraction
    library("BSgenome.Hsapiens.UCSC.hg19")
    tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep = ":"))
    genome <- eval(parse(text = tmp))

    seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
    sequences <- lapply(positionsList, function(x){ seqExtract(positions = x, genome = genome, seq_range) })

    #Write splitted sequences
    tempDir <- paste(outname, "_temp", sep = "")
    dir.create(tempDir)

    seqs <- sequences$target
    target_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "target")
    seqs <- sequences$random
    random_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "random")
    rm(sequences)
    rm(seqs)
    invisible(replicate(3, gc()))

    #motif search
    target_positionsList <- splitSeqMotDist(filenames = target_all_filenames, motif_list = motif, min.score = min.score)
    file.remove(target_all_filenames)
    random_positionsList <- splitSeqMotDist(filenames = random_all_filenames, motif_list = motif, min.score = min.score)
    file.remove(random_all_filenames)

    #computation of enrichment score
    target_mot_posi <- unlist(target_positionsList)
    ctrl_mot_posi <- unlist(random_positionsList)

    #motif enrichment score
    enrichment_scores <- enrichScoreDist(
        target_mot_posi = target_mot_posi,
        ctrl_mot_posi = ctrl_mot_posi,
        seq_range = seq_range,
        motif_name = motif_names,
        nDMP_IDs = nDMP_IDs,
        outname = outname,
        plot_draw = FALSE
    )
    All_enrichments <- c(All_enrichments, list(enrichment_scores))

    file.remove(tempDir)
}

#Vidualization
ranks <- seq(seq_range[1], seq_range[2], length = 1001) # x-axis values
data <- sapply(All_enrichments, cbind)
rownames(data) <- ranks
colnames(data) <- comp_names

data.df <- as.data.frame(data)
data.df %>%
    mutate(Position = rownames(data.df)) %>%
    mutate_at(vars("Position"), as.integer) -> data.df

data.df2 <- tidyr::pivot_longer(data.df, cols = -Position, names_to = "comparison", values_to = "ES")

##ggplot
data_names <- names(uninherit_demet_list)

cols <- brewer.pal(6, "Dark2")
g_list <- NULL
for(i in seq(length(data_names))){
    data.df3 <- filter(data.df2, comparison == eval(data_names[i]))

    g <- ggplot(data.df3, aes(x = Position, y = ES, color = eval(cols[i])))
    g <- g + geom_line()
    g <- g + ylim(-0.1, 0.25)
    g <- g + xlab("Distance from CpG") # X-axis label
    g <- g + ylab("Enrichment Score") # Y-axis label
    g <- g + geom_hline(yintercept = 0) # wtire a horizontal line at intercept=0
    g <- g + geom_vline(xintercept = 0, color = "Gray51") # wtire a horizontal line at intercept=0
    g <- g + scale_colour_manual(values = c(cols[i]))
    ## Settings of axes
    g <- g + theme(
        axis.text.x = element_text(color = "black", size = 14), # X-axis texit setting
        axis.text.y = element_text(color = "black", size = 14), # Y-axis texit setting
        axis.ticks = element_line(color = "black", size = .5), # Setting of scale ticks line
        axis.title = element_text(color = "black", size = 14), # Setting of scale ticks texit
        panel.background = element_rect(fill = "white", color = NA), # panel background color
        panel.border = element_rect(fill = NA, color = "black"), # panel border
        panel.grid.major = element_blank(), # panel major grid
        panel.grid.minor = element_blank() # panel minor grid
        )
    ## Settings of legend
    g <- g + theme(legend.position = 'none')
    g_list <- c(g_list, list(g))
}

g_all <- plot_grid(plotlist = g_list, nrow = 1)
ggsave(g_all, file = "uninherited_enrichment_plot.pdf", width = 35, height = 7)

#===================================================================================================

#  Enrichment heatmap
#
# 1. GATA6 ChIPmeantation (Fig. 4E)
# 2. ATAC-seq (Fig. 5B)

#===================================================================================================
#-------------------------------------------------
#bed file expoert
#-------------------------------------------------
extension <- c(0, 1)
ui.1_bed_list <- NULL
for(i in 2:6){
    ui.1_bed_list  <- c(
        ui.1_bed_list ,
        list(dmr2bed(
            data_matrix = selDataMatrix,
            ControlColnum = 1,
            TreatmentColnum = i,
            extension = extension
        ))
    )
}

names(ui.1_bed_list) <- c("Demet.0h_48h", "Demet.0h_54h", "Demet.0h_60h", "Demet.0h_66h", "Demet.0h_72h")


dmr.probeID <- lapply(ui.1_bed_list, function(x){ x[, 4] })


#uninherited
uninherit_demet_list <- list(
    Demet.0h_48h = ui.1_bed_list[[1]],
    Demet.0h_54h = ui.1_bed_list[[2]][!(dmr.probeID[[2]] %in% unlist(dmr.probeID[1])),],
    Demet.0h_60h = ui.1_bed_list[[3]][!(dmr.probeID[[3]] %in% unlist(dmr.probeID[1:2])),],
    Demet.0h_66h = ui.1_bed_list[[4]][!(dmr.probeID[[4]] %in% unlist(dmr.probeID[1:3])),],
    Demet.0h_72h = ui.1_bed_list[[5]][!(dmr.probeID[[5]] %in% unlist(dmr.probeID[1:4])),]
)


#export uninherited demet beds
dir.create("data/Methyl850_DE/bed_d1", recursive = TRUE)

for(j in 1:5){ write.table(
    x = uninherit_demet_list[[j]],
    file = paste0("data/Methyl850_DE/bed_d1/", names(uninherit_demet_list)[j], ".bed"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)}

#-------------------------------------------------
#heatmap of GATA6 ChIPmentation
#-------------------------------------------------
bed_path <- "data/ethyl850_DE/d1/"
beds <- c("Demet.0h_48h.bed", "Demet.0h_54h.bed", "Demet.0h_60h.bed", "Demet.0h_66h.bed", "Demet.0h_72h.bed")

#bigwig files
bws_ch <- c("GATA6_0h.sorted.bw", "GATA6_48h.sorted.bw", "GATA6_54h.sorted.bw", "GATA6_60h.sorted.bw", "GATA6_66h.sorted.bw", "GATA6_72h.sorted.bw")
bws_ac <- c("atac_merged_0h.bw", "atac_merged_48h.bw", "atac_merged_54h.bw", "atac_merged_60h.bw", "atac_merged_66h.bw", "atac_merged_72h.bw")

titles <- c("0 h", "48 h", "54 h", "60 h", "66 h", "72 h")

hm_list_all_ChIP <- list()
hm_list_all_ATAC <- list()
for(k in seq(length(beds))){
    #ChIPmenatation
    hm_list_ChIP <- NULL
    for(l in seq(lHeatmapength(bws_ch))){
        hm_chip <- covHeatmap(
            target = paste0(bed_path, beds[k]),
            bw = paste0(bw_path, bws_ch[l]),
            add_profile = FALSE,
            heatmap_legend_param = list(title = titles[l], title_gp = gpar(fontsize = 12, fontfamily = "HersheySans")),
            column_title = titles[l],
            column_title_gp = gpar(fontsize = 12, fontfamily = "HersheySans")
            col = colorRamp2(c(0, 15), c("darkblue", "darkgoldenrod1"))
        )
        hm_list_ChIP <- c(hm_list_ChIP, list(hm_chip))
    }
    #integrate heatmaps of a same DMRs
    hm_list_ChIP_int <- hm_list_ChIP[[1]] + 
        hm_list_ChIP[[2]] + 
        hm_list_ChIP[[3]] + 
        hm_list_ChIP[[4]] + 
        hm_list_ChIP[[5]] + 
        hm_list_ChIP[[6]]
    hm_list_all_ChIP <- c(hm_list_all_ChIP, list(hm_list_ChIP_int))
    
    #ATAC-seq
    hm_list_ATAC <- NULL
    for(m in seq(length(bws_ac))){
        hm_atac <- covHeatmap(
            target = paste0(bed_path, beds[k]),
            bw = paste0(bw_path, bws_ac[m]),
            add_profile = FALSE,
            heatmap_legend_param = list(title = titles[m], title_gp = gpar(fontsize = 12, fontfamily = "HersheySans")),
            column_title = titles[m],
            column_title_gp = gpar(fontsize = 12, fontfamily = "HersheySans"),
            cols = c("white", "red")
        )
        hm_list_ATAC <- c(hm_list_ATAC, list(hm_atac))
    }
    #integrate heatmaps of a same DMRs
    hm_list_ATAC_int <- hm_list_ATAC[[1]] + 
        hm_list_ATAC[[2]] + 
        hm_list_ATAC[[3]] + 
        hm_list_ATAC[[4]] + 
        hm_list_ATAC[[5]] + 
        hm_list_ATAC[[6]]
    hm_list_all_ATAC <- c(hm_list_all_ATAC, list(hm_list_ATAC_int))
}

# integrated heatmaps of different DMRs
grab_list_ChIP <- lapply(hm_list_all_ChIP, function(x){grid.grabExpr(draw(x))})
grab_list_AtAC <- lapply(hm_list_all_ATAC, function(x){grid.grabExpr(draw(x))})

#plot
g_C <- plot_grid(plotlist = grab_list_ChIP, ncol = 1, align = "hv", axis =  "tblr")
g_ATAC <- plot_grid(plotlist = grab_list_ATAC, ncol = 1, align = "hv", axis =  "tblr")

#===================================================================================================

#expression analysis for demet probes (Fig. 5A)

#===================================================================================================
#bed files of demethylated region
#used for ctss count (the ctss count data are stored in data/Methyl850_DE/exp_table)
extension <- c(-250, 250)
ui.250_bed_list <- NULL
for(i in 2:6){
    ui.250_bed <- dmr2bed(data_matrix = selDataMatrix, ControlColnum = 1, TreatmentColnum = i, extension = extension)
    ui.250_bed_list <- c(ui.250_bed_list, list(ui.250_bed))
}
names(ui.250_bed_list) <- c("DE_0.48", "DE_0.54", "DE_0.60", "DE_0.66", "DE_0.72")

dir.create("data/Methyl850_DE/bed_d500", recursive = TRUE)
for(i in seq(ncol(ui.250_bed_list))){
    write.table(
        ui.250_bed_list[[j]],
        file = paste0("data/Methyl850_DE/bed/", names(ui.250_bed_list)[j], ".bed"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
    )
}

#read the xpression tabels
data_dir <- "data/Methyl850_DE/exp_table"
files <- c(
    "iPS-DE_0.48_demethyl.TPM.txt",
    "iPS-DE_0.54_demethyl.TPM.txt",
    "iPS-DE_0.60_demethyl.TPM.txt",
    "iPS-DE_0.66_demethyl.TPM.txt",
    "iPS-DE_0.72_demethyl.TPM.txt"
)
sample_names <- sapply(files, function(x){ gsub(".TPM.txt", "", x) })

#expression table list
exp_list <- NULL
list_names <- NULL
for(i in seq(length(files))){
    exp_data <- read.table(
        file = paste0("data/Methyl850_DE/exp_table/", files[i]),
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    #column reorder
    id_index <- 1
    reorder_index <- order(colnames(exp_data)[-id_index])
    exp_data2 <- cbind(probeID = exp_data[, id_index], exp_data[, - id_index][, reorder_index])

    #output list
    exp_list <- c(exp_list, list(exp_data2))
}
names(exp_list) <- sample_names

#mean and sd of rxpression
group_ids <- c("iPS_DE_00h","iPS_DE_48h","iPS_DE_54h","iPS_DE_60h","iPS_DE_66h","iPS_DE_72h")

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


    colnames(mean_expMatrix) <- paste0(group_ids, "exp_mean")
    rownames(mean_expMatrix) <- exp_table$probeID
    mean_expMatrix_list <- c(mean_expMatrix_list, list(mean_expMatrix))

    colnames(sd_expMatrix) <- paste0(group_ids, "exp_sd")
    rownames(sd_expMatrix) <- exp_table$probeID
    sd_expMatrix_list <- c(sd_expMatrix_list, list(sd_expMatrix))
}

names(mean_expMatrix_list) <- sample_names
names(sd_expMatrix_list) <- sample_names

#Vidualization
g_list <- NULL
for(i in seq(length(mean_expMatrix_list))){
    df <- as.data.frame(mean_expMatrix_list[[i]])
    names(df) <- c(0, 48, 54, 60, 66, 72)

    g <- tcLinePlot(df = df, is.log = FALSE, color = "green4", fill = alpha("palegreen1", 0.4), y_lim = c(-9, 15))
    glist <- c(g_list, list(g))
}

g_all <- plot_grid(plotlist = g_list, nrow = 1)
ggsave(g_all, file = "mean.tpm_heatmap_demet.pdf", width = 42, height = 7)

#===================================================================================================

#Functional Annotation of DMRs (Fig. S6A)

#===================================================================================================
#Processing of gencode gff3 file
filepath <- "data/Methyl850_DE/gencode.v19.gene_only.gff3"
gencode.gff <- readGFF(filepath, version = 1)
gencode.gr <- makeGRangesFromDataFrame(gencode.gff, keep.extra.columns = TRUE)

annotations <- strsplit(as.character(gencode.gr$group), ";")

gene_name <- sapply(annotations, function(x){
    gsub("gene_name=", "", x[6])
})

gencode.gr$gene_name <- gene_name

#Promoter extraction
upstream <- 1000
downstream <- 200
prom.gr <- promoters(gencode.gr, upstream = upstream, downstream = downstream)

#loading  enhances
eh.gr <- import("data/human_permissive_enhancers_phase_1_and_2.bed")

#count overlaps
Promoter <- NULL
Enhancer <- NULL
Both <- NULL
Unannotated <- NULL
for(i in seq(length(uninherit_demet_list))){
    prom_count <- countOverlaps(uninherit_demet_list[[i]], prom.gr)
    eh_count <- countOverlaps(uninherit_demet_list[[i]], eh.gr)

    nboth <- sum(prom_count != 0 & eh_count != 0)
    weight_both <- nboth / length(uninherit_demet_list[[i]])
    Both <- c(Both, weight_both)

    nprom <- sum(prom_count != 0) - nboth
    weight_prom <- nprom / length(uninherit_demet_list[[i]])
    Promoter <- c(Promoter, weight_prom)

    neh <- sum(eh_count != 0) - nboth
    weight_eh <- neh / length(uninherit_demet_list[[i]])
    Enhancer <- c(Enhancer, weight_eh)

    nunannotated <- length(uninherit_demet_list[[i]]) - nboth - nprom - neh
    weight_nunannotated <- 1 - weight_both - weight_prom - weight_eh
    Unannotated <- c(Unannotated, weight_nunannotated)
}
df <- data.frame(comp = names(uninherit_demet_list), Promoter, Enhancer, Unannotated)

#vidualization
theme <- theme(panel.background = element_blank(), # initialization
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.text.x = element_text(size = 10, angle = 90),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10),
            axis.ticks = element_line(size = 0.5),
            axis.line = element_line(size = 0.5),
            plot.margin = unit(c(1, 1, 1, 1), "line")
    )


df2 <- pivot_longer(df, col = -comp, names_to = "Annotation", values_to = "Proportion")
g <- ggplot(df2, aes(x = comp, y = Proportion, fill = Annotation))
g <- g + geom_bar(stat = "identity")
g <- g + scale_fill_manual(values = c("royalblue4", "royalblue1", "lightskyblue1"))
g <- g + theme
g

outprefix <- "annotation_proportion.pdf"
ggsave(plot = g, file = outprefix, width = 7, height = 7)
