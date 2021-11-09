#===================================================================================================#

#Preprocessing

#===================================================================================================#
selDataMatrix <- lumiMethyNorm(
    fileName = "Matrix_signal_intensities.txt",
    inputtype = "signal",
    sample_names = c("Hep_D0", "Hep_D7", "Hep_D14", "Hep_D21", "Hep_D28")
)

#===================================================================================================

# correlation matrix (Fig. 1C)

#===================================================================================================
cor.exp <- cor(selDataMatrix)

hc <- hclust(dist(cor.exp))
breaks1 <- seq(0.95, 1, length = 500) #the range should be optimized
mycols1 <- colorRamp2(breaks1, colorRampPalette(brewer.pal(9, "GnBu"))(500)) # color for motif enrichment

ht1 <- Heatmap(cor.exp,
                col = mycols1,
                name = "Correlation Coefficient", #title of legend, 
                row_names_gp = gpar(fontsize = 12), # Text size for row names
                column_names_gp = gpar(fontsize = 12), # Text size for column names
                column_names_rot = 90,
                heatmap_legend_param = list(title = "R",
                    title_gp = gpar(fontsize = 12),
                    labels_gp = gpar(fontsize = 6),
                    grid_height = unit(0.8, "cm"),
                    grid_width = unit(0.8, "cm"),
                    title_position = c("topcenter")
                    ),
                row_dend_width = unit(2, "cm"), #row dendrogram width
                column_dend_height = unit(2, "cm"), #column dendrogram height
                row_dend_side = "right",
                column_dend_side = "bottom",
                cluster_rows = hc,
                cluster_columns = hc,
                cell_fun = function(j, i, x, y, width, height, fill){
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
                }
            )
ht1
pdf("Hep_correlationplot.pdf")
draw(ht1)
dev.off()

#===================================================================================================

# DMR b2b plot (Fig. 1D)

#===================================================================================================
#Count DMP
comp <- data.frame(D0_D7=c(1,2), D7_D14=c(2,3), D14_D21=c(3,4), D21_D28=c(4,5)) # comparisons

M_ndiffprobes <- sapply(comp, function(x){
    sum(selDataMatrix[,as.numeric(x[1])] - selDataMatrix[,x[2]] <= -2)
})
D_ndiffprobes <- sapply(comp, function(x){
    sum(selDataMatrix[,as.numeric(x[1])] - selDataMatrix[,x[2]] >= 2)
})

nprobes_df <- data.frame(Time=c("D0-D7", "D7-D14", "D14-D21", "D21-D28"), Methylated=M_ndiffprobes, Demethylated=D_ndiffprobes)

#plot 
nprobes_b2b_df <- tidyr::pivot_longer(nprobes_df, col = -Time, , names_to = "MetDemet", values_to = "nprobes")


the_order <- rev(as.character(unlist(nprobes_b2b_df %>%
    filter (MetDemet == "Methylated") %>%
    select(Time)
)))

theme <- theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(size = 10,hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1),"line"),
    panel.grid.major.x = element_line(linetype = "dotted", size = 0.3, color = "#3A3F4A")
)

g <- ggplot(nprobes_b2b_df, aes(x=Time))
g <- g + scale_x_discrete(limits = the_order)

#Methylated plot (left)
g.M <- g + geom_bar(
    data=nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]]=="Methylated",],
    aes(y=nprobes), fill=rgb(0.8, 0, 0.2, alpha=0.9),
    stat = "identity"
)
g.M <- g.M + scale_y_reverse(limits=c(max(nprobes_b2b_df[["nprobes"]])*1.2,0))
g.M <- g.M + ggtitle('Methylated')
g.M <- g.M + ylab("Number of probes")
g.M <- g.M + coord_flip()
g.M <- g.M + theme
g.M <- g.M + theme(
    axis.text.x = element_text(size = 10), 
    axis.ticks.x = element_line(size=0.5), 
    axis.line.x = element_line(size=0.5)
)

#Demethylated plot (Right)
g.D <- g + geom_bar(
    data=nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]]=="Demethylated",],
    aes(y=nprobes), fill=rgb(0, 0.7, 0, alpha=0.9), 
    stat = "identity"
)
g.D <- g.D + scale_y_continuous(limits=c(0,max(nprobes_b2b_df[["nprobes"]])*1.2))
g.D <- g.D + ggtitle('Demethylated')
g.D <- g.D + ylab("Number of probes")
g.D <- g.D + coord_flip()
g.D <- g.D + theme
g.D <- g.D + theme(
    axis.text.x = element_text(size = 10), 
    axis.ticks.x = element_line(size=0.5), 
    axis.line.x = element_line(size=0.5)
)


#labels (Center)
g.time <- g + geom_bar(
    data=nprobes_b2b_df[nprobes_b2b_df[["MetDemet"]]=="Demethylated",],
    aes(y=0, fill = 'white'), 
    alpha=0, stat = "identity"
)
g.time <- g.time + geom_text( aes( y = 0,  label = as.character(Time)), size = 3)
g.time <- g.time + ggtitle("Time")
g.time <- g.time + coord_flip()
g.time <- g.time + ylab("")
g.time <- g.time + theme
g.time <- g.time + theme(panel.grid.major.x = element_blank()) 

g_all <- grid.arrange(
    g.M, g.time, g.D,
    ncol=3, 
    widths=c(4.25/10,1.5/10,4.25/10)
)

ggsave(g_all, file="nDMP_b2bplot.pdf", width=7, height=7)

#===================================================================================================

# Functional analysis of DMRs (Fig. S2)

    #GO analysis using eGREAT and then summaryze using rrvgo

#===================================================================================================
COMP <- data.frame(Hep_00_07 = c(1, 2), Hep_07_14 = c(2, 3), Hep_14_21 = c(3, 4), Hep_21_28 = c(4, 5))
topn <- 5

for(i in c("Methyl, Demethyl")){
    #Construction of Entrez ID list
    enrich_reducedTerms_BP.df <- NULL
    for(j in seq(length(COMP))){
        job <- DMPrGREAT(
            data = selDataMatrix,
            out_prefix = paste0(i, "_", names(COMP)[j]),
            ControlColnum = COMP[1, j],
            TreatmentColnum = COMP[2, j],
            MethylDemethyl = i,
        )

        #extraction of enriched ensemble gene table
        tb = getEnrichmentTables(job)

        #Calculation of the similarity matrix
        simMatrix <- calculateSimMatrix(
            tb$"GO Biological Process"$ID,
            orgdb = "org.Hs.eg.db",
            ont = "BP",
            method = "Rel"
        )
        scores <- setNames(-log10(tb$"GO Biological Process"$Hyper_Adjp_BH), tb$"GO Biological Process"$ID)
        reducedTerms <- reduceSimMatrix(
            simMatrix = simMatrix,
            scores = scores,
            threshold = threshold,
            orgdb = "org.Hs.eg.db"
        )
        reducedTerms2 <- left_join(reducedTerms, tb$"GO Biological Process", by = c("go" = "ID"))

        #summarize by cluster
        ndmp <- length(DmpId(selDataMatrix, COMP[1, j], COMP[2, j], MethylDemethyl = i))

        reducedTerms2 %>%
            group_by(parentTerm) %>%
            summarise(total_Region_Hits = sum(Hyper_Foreground_Region_Hits), total_score = sum(score)) %>%
            mutate(Hit_ratio = total_Region_Hits / ndmp) %>%
            mutate(sample = names(COMP)[j]) -> parents_summry

        parents_summry_top <- parents_summry[order(parents_summry$total_score, decreasing = TRUE)[seq(topn)],]
        enrich_reducedTerms_BP.df <- rbind(enrich_reducedTerms_BP.df, parents_summry_top)
    }

    enrich_reducedTerms_BP.df$Reduced_GO <- factor(
        enrich_reducedTerms_BP.df$parentTerm,
        levels = rev(unique(enrich_reducedTerms_BP.df$parentTerm))
    )

    #plot
    g <- ggplot(data = enrich_reducedTerms_BP.df, aes(x = sample, y = Reduced_GO, color = total_score, size = Hit_ratio))
    g <- g + geom_point()
    g <- g + scale_color_gradient(low = "blue", high = "red")
    g <- g + theme_light()
    g <- g + theme(axis.text.x = element_text(size = 15, angle = 90),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
        )
    ggsave(g, file = paste0("rGREATrrvgo_", i, "__GO.BP_summary.pdf"), width = 14, height = 10.5)
}

#===================================================================================================

# Heatmap for methylation regulating TF (Fig. 2B, D)

#===================================================================================================
#Motif overrepresentation analysis
motifDBList <- IMAGE_PWMlist
comp <- c("Hep_00_07", "Hep_07_14", "Hep_14_21", "Hep_21_28")
    
for(i in c("Methyl", "Demethyl")){
    for(j in 1:4){
        outdir <- paste0("data/Methyl850_HEP/motif_analysis/", comp[j])
        dir.create(outdir, recursive = TRUE)
        
        MotScr(
            infile = selDataMatrix,
            motifDBList = motifDBList,
            outname = paste0(outdir, "/", ),
            ControlColnum = j,
            TreatmentColnum = j + 1,
            MethylDemethyl = i
            )
    }
}

#Heatmap
for(i in c("Methyl", "Demethyl")){
    #load the data
    data_list <- lapply(comp, function(x){
        read.table(paste0("data/Methyl850_HEP/motif_analysis/", x, "/", i, "_mot_analysis_result.txt"),
        sep = "\t",
        header = T,
        stringsAsFactors = FALSE)
    })
    names(data_list) <- comp
    
    #Processing of motif names
    motif_names <- sapply(sapply(data_list[[1]][, 1], function(x){ strsplit(as.character(x), "\\.") }), function(x){ x[1] })
    
    #P-value matrix
    pvals <- as.data.frame(lapply(data_list, function(x){ x[, 7] }))
    rownames(pvals) <- motif_names
    colnames(pvals) <- paste0(comp, "_p")
    
    #Conversion of pvalue = 0 to 1.00e-300
    pvals.matrix <- as.matrix(pvals)
    pvals.matrix[pvals.matrix == 0] <- 1.00e-300
    pvals <- as.data.frame(pvals.matrix)
    
    # "significant" matrix
    sig_matrix <- sapply(data_list, function(x){ x[, 15] })
    rownames(sig_matrix) <- motif_names
    colnames(sig_matrix) <- paste0(comp, "_sig")
    
    #Conbine "Significant" matrix and p-value matrix as.data frame
    pval_sig <- data.frame(motif = rownames(pvals), as.data.frame(pvals), as.data.frame(sig_matrix))
    
    #Selection of rows which has at least one pvalue <= 10e-50 and at least one "significant"
    pcutoff <- 1
    
    pval_sig %>%
        select(ends_with("_p")) %>%
        apply(., 1, min) < pcutoff -> pval_index
    pval_sig %>%
        select(ends_with("_sig")) %>%
        apply(., 1, function(x){ any(x == "significant") }) -> sig_index
    
    sel_pval_sig <- filter(pval_sig, pval_index & sig_index)
    
    #non-significant p-value to NA
    sel_pval_sig %>%
        select(ends_with("_p")) -> sel_pvals
    sel_pval_sig %>%
        select(ends_with("_sig")) ->  sel_sigs
    sel_pvals[sel_sigs != "significant" | is.na(sel_sigs)] <- "NA"
    
    sel_pval_sig_na <- cbind(sel_pvals, sel_sigs)
    
    #load CAGE expression table
    load("data/CAGE_HEP/gene_expMatrix.RData")
    
    E_TC_names <- c("Hep_00", "Hep_07", "Hep_14", "Hep_21", "Hep_28")
    colnames(mean.gene_expMatrix) <- paste0(E_TC_names, "_exp")
    colnames(sd.gene_expMatrix) <- paste0(E_TC_names, "_sd")
    mean_sd.gene_expMatrix <- cbind(mean.gene_expMatrix, sd.gene_expMatrix)
    
    #integrate motif enrichment and CAGE expression data
    gene_names <- as.vector(t(as.data.frame(strsplit(as.character(rownames(sel_pval_sig_na)), "_")))[, 1])
    
    #selection of mean.gene.expMatrix
    sel_mean_sd.gene_expMatrix <- t(sapply(gene_names, function(x){
        if(any(rownames(mean_sd.gene_expMatrix) %in% x)){
            mean_sd.gene_expMatrix[rownames(mean_sd.gene_expMatrix) %in% x,]
        } else{
        return(rep(NA, ncol(mean_sd.gene_expMatrix)))
        }
    }))
    
    sel_all_df <- cbind(sel_pval_sig_na, sel_mean_sd.gene_expMatrix)
    sel_all_df <- mutate(sel_all_df, motif = rownames(sel_all_df))
    
    #remove low expression motifs
    TPMcutoff <- 50
    
    highExp_sel_all_df <- NULL
    for(i in seq(length(comp))){
        Msig_index <- sel_all_df[, paste0(comp[i], "_sig")] == "significant"
        TPM_index <- as.vector(apply(sel_all_df[, paste0(E_TC_names[c(i, i + 1)], "_exp")], 1, max) >= TPMcutoff)
        highExp_sel_all_df <- rbind(highExp_sel_all_df, sel_all_df[which(Msig_index & TPM_index),])
    }
    highExp_sel_all_df <- highExp_sel_all_df %>% distinct(motif, .keep_all = TRUE)
    
    #plot heatmap
    #Extract pvals for heatmap
    highExp_sel_pvals <- select(highExp_sel_all_df, ends_with("_p"))
    rownames(highExp_sel_pvals) <- highExp_sel_all_df[, "motif"]
    colnames(highExp_sel_pvals) <- gsub("_p", "", colnames(highExp_sel_pvals))
    
    #log10 conversion
    log_highExp_sel_pvals <- apply(highExp_sel_pvals, c(1, 2), function(x){-log(as.numeric(x), 10) })
    
    #heatmap
    df1 <- log_highExp_sel_pvals
    breaks1 <- seq(0, 150, length = 500)
    mycols1 <- colorRamp2(breaks1, colorRampPalette(brewer.pal(9, "GnBu"))(500)) # color for motif enrichment
    hc.rows <- hclust(dist(as.matrix(replace(df1, is.na(df1), 0))))
    
    ## cluster color
    if(i == "Methyl"){
        clust_num <- 4  
    } else{
        clust_num <- 5 
    }
    
    clusters <- cutree(hc.rows, clust_num) #cluster annotation
    clust_col <- brewer.pal(clust_num, "Pastel1") #cluster colors
    names(clust_col) <- seq(clust_num)
    ha_row = rowAnnotation(df = data.frame(cluster = clusters), col = list(cluster = clust_col), width = unit(1, "cm"))
    
    ht1 <- Heatmap(df1,
        width = unit(6, "cm"), #With of heatmap
        name = "-10g10P", #title of legend
        column_title = "diff probe nubers and Motif enrichment",
        column_title_gp = gpar(fontsize = 10),
        row_title = NA,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 5), # Text size for row names
        column_names_gp = gpar(fontsize = 8), # Text size for column names
        col = mycols1,
        cluster_rows = hc.rows,
        row_dend_width = unit(4, "cm"), #row dendrogram width
        row_dend_reorder = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lty = 1, lwd = 0.5),
        na_col = "grey50",
        cell_fun = function(j, i, x, y, width, height, fill){
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
        },
        right_annotation = ha_row
    )
    
    pdf(paste0(i , "_motif_enirch_heatmap.pdf"))
    draw(ht1)
    dev.off()
    
    save(list = ls(), file = paste0("data/Methyl_850_HEP/", i, "_motif_heatmap.RData"))
}
#===================================================================================================

# exprssion line plot (Fig. 2C and 2E; Fig. S3)

#===================================================================================================
for(i in c("Methyl", "Demethyl")){
    #load tmp data
    load(paste0("data/Methyl850_HEP/", i ,"_motif_heatmap.RData"))
    
    ## add acluster metadata
    highExp_sel_all_df2 <- cbind(highExp_sel_all_df, cluster = clusters)
    
    # motif cluster expression and sd
    highExp_sel_exp_sd <- select(highExp_sel_all_df2, motif, cluster, ends_with("_exp"), ends_with("_sd"))
    
    if(i == "Methyl"){
        cluster_bum <- 4
    } else{
        cluster_num <- 5
    }
    
    y.df_list <- NULL
    for(j in seq(clust_num)){
        highExp_sel_exp_sd %>%
            filter(cluster == j) %>%
            arrange(motif) -> cluster_gene_exp
            
        #extract gene names
        gene <- sapply(strsplit(as.character(unlist(select(cluster_gene_exp, motif))), "_"), function(x){ x[1] })
        
        #remove dupicated genes
        cluster_gene_exp2 <- cluster_gene_exp[which(!duplicated(gene)),]
        gene2 <- gene[which(!duplicated(gene))]
        
        #extract only expression data and formatting for plot
        df2 <- as.data.frame(cbind(Gene = gene2, select(cluster_gene_exp2, ends_with("_exp"))))
        colnames(df2)[-1] <- c(0, 7, 14, 21, 28)
        
        y <- tidyr::pivot_longer(df2, cols = -Gene, names_to = "Day", values_to = "TPM")
        
        #extract only sd data and formatting for plot
        df2.sd <- as.data.frame(cbind(Gene = gene2, select(cluster_gene_exp2, ends_with("_sd"))))
        colnames(df2.sd)[-1] <- c(0, 7, 14, 21, 28)
        sd.y <- tidyr::pivot_longer(df2.sd, cols = -Gene, names_to = "Day", values_to = "SD")
        
        #integrate TPM and SD
        y.df <- cbind(y, select(sd.y, SD))
        y.df_list <- c(y.df_list, list(y.df))
    }
    
        #plot
    theme <- theme(panel.background = element_blank(), # initialization
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(colour = "black", size = 9),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(colour = "black", size = 9),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "bottom",
        plot.margin = unit(c(1, 1, 1, 1), "line")
    )
    
    g_list <- NULL
    for(k in seq(length(y.df_list))){
        df3 <- y.df_list[[k]]
        cluster_name <- paste0("Cluster ", k)
        
        ## Color setting
        max.tpm <- df3 %>% group_by(Gene) %>% summarise(max = max(TPM))
        max.tpm <- floor(max.tpm[, "max"])
        TPM_ranged_col <- colorRampPalette(brewer.pal(9, "YlGn"))(max(max.tpm))
        mycols <- sapply(max.tpm, function(x){ TPM_ranged_col[x] })
        
        #For legend of max TPM (color of each gene)
        dat = data.frame(max_TPM = seq(max(max.tpm)), y = 1, color = TPM_ranged_col)
        gc <- ggplot(dat, aes(max_TPM, y, colour = max_TPM)) +
                geom_point(size = 0) +
                scale_colour_gradientn(colors = TPM_ranged_col) +
                theme(legend.text = element_text(colour = "black", size = 9))
        
        tpm_color_legend <- get_legend(gc)
        
        ## ggplot2
        g <- ggplot()
        g <- g + geom_line(df3, mapping = aes(x = as.numeric(Day), y = TPM, color = Gene, group = Gene), size = 0.6) # draw lines
        g <- g + geom_point(df3, mapping = aes(x = as.numeric(Day), y = TPM, color = Gene, group = Gene), size = 0.9) #draw points
        g <- g + geom_errorbar(df3, mapping = aes(x = as.numeric(Day), ymin = TPM - SD, ymax = TPM + SD, color = Gene, group = Gene), width = 0.2, size = 0.3) #error bar
        g <- g + ggtitle(cluster_name) #title
        g <- g + xlab("Day") #y axis label
        g <- g + scale_color_manual(values = as.vector(mycols)) #color of each data
        g <- g + scale_x_continuous(breaks = c(0, 7, 14, 21, 28))
        g <- g + theme_bw(10) #basic theme
        g <- g + theme #custam theme
        g <- g + guides(color = guide_legend("Gene"), fill = guide_legend("Gene")) #lege
        
        g_legend <- get_legend(g)
        
        g2 <- plot_grid(g + theme(legend.position="none"), tpm_color_legend, nrow=1, ncol=2, rel_widths = c(9, 1))
        g3 <- plot_grid(g2, g_legend, nrow=2, ncol=1, rel_heights = c(7.5, 2.5)) 
        
        g_list <- c(g_list, list(g3))
    }
    
    g4 <- plot_grid(plotlist = g_list, nrow = 2, align = "hv")
    
    ggsave(g4, file = paste0(i, "_TF_exp.pdf"), width = 21, height = 14)
}
