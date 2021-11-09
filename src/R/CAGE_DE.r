#===================================================================================================#

# Marker gene exprssion (Fig. S5A)
#
#
#   description: This script vidualize expression of indicated genes as a line plot.
#===================================================================================================#
#loding the CAGE expression data
expMatrix <- read.table("iPS_DE_diff_F5.TPM.txt", header = TRUE, row.name = NULL, sep = "\t", stringsAsFactors = FALSE)
expMatrix <- expMatrix[, c(1, 3, 10, 17, 5, 12, 19, 8, 15, 22, 7, 14, 21, 2, 9, 16, 6, 13, 20)] # reorder based on time point

#Compute mean and sd of replicates
mean.expMatrix <- cbind(
    apply(expMatrix[, c(2:4)], 1, mean),
    apply(expMatrix[, c(5:7)], 1, mean),
    apply(expMatrix[, c(8:10)], 1, mean),
    apply(expMatrix[, c(11:13)], 1, mean),
    apply(expMatrix[, c(14:16)], 1, mean),
    apply(expMatrix[, c(17:19)], 1, mean)
)
colnames(mean.expMatrix) <- c(0, 48, 54, 60, 66, 72)
rownames(mean.expMatrix) <- expMatrix$prmtrID

##compute sd of replicates
sd.expMatrix <- cbind(
    apply(expMatrix[, c(2:4)], 1, sd),
    apply(expMatrix[, c(5:7)], 1, sd),
    apply(expMatrix[, c(8:10)], 1, sd),
    apply(expMatrix[, c(11:13)], 1, sd),
    apply(expMatrix[, c(14:16)], 1, sd),
    apply(expMatrix[, c(17:19)], 1, sd)
)
colnames(sd.expMatrix) <- c(0, 48, 54, 60, 66, 72)
rownames(sd.expMatrix) <- expMatrix$prmtrID

#extraction of GOI
psc_marker <- c("POU5F1", "NANOG")
ps_marker <- "T"
de_marker <- c("SOX17",  "FOXA2", "GATA6")

gois <- c(psc_marker, ps_marker,de_marker)
for(i in gois){
    goi <- i

    gene_names <- sapply(strsplit(rownames(mean.expMatrix), "@"), function(x){x[2]})

    goi.mean.expMatrix <- mean.expMatrix[gene_names %in% goi,, drop = FALSE]
    goi.sd.expMatrix <- sd.expMatrix[gene_names %in% goi,, drop = FALSE]

    ##order based on promoter number
    goi.mean.expMatrix <- goi.mean.expMatrix[order(as.numeric(gsub("p", "", sapply(strsplit(rownames(goi.mean.expMatrix), "@"), function(x){ x[1] })))),]
    goi.sd.expMatrix <- goi.sd.expMatrix[order(as.numeric(gsub("p", "", sapply(strsplit(rownames(goi.sd.expMatrix), "@"), function(x){ x[1] })))),]

    #remove low expreeesion promoters
    cutoff <- 5
    tpm_cutoff <- which((apply(goi.mean.expMatrix, 1, max) >= cutoff))
    goi.mean.expMatrix <- goi.mean.expMatrix[tpm_cutoff,, drop = FALSE]
    goi.sd.expMatrix <- goi.sd.expMatrix[tpm_cutoff,, drop = FALSE]

    if(nrow(goi.mean.expMatrix) == 0) next

    library(ggplot2)
    library(tidyverse)
    library(reshape2)
    library(ggsci)

    goi.mean.expMatrix %>%
        as.data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        pivot_longer(cols = -Gene, values_to = "TPM", names_to = "Hours") -> y

    goi.sd.expMatrix %>%
        as.data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        pivot_longer(cols = -Gene, values_to = "sd", names_to = "Hours") -> sd.y

    y.df <- cbind(y, sd = sd.y[, 3])
    y.df$Hours <- as.numeric(y.df$Hours, levels = c(0, 48, 54, 60, 66,72))

    theme <- theme(
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(color = "black", size = 20),    # axis texit setting
        axis.ticks = element_line(color = "black", size = .5),    # Setting of scale ticks line
        axis.title = element_text(color = "black", size = 20),    # Setting of scale ticks texit
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=20),
        plot.margin = unit(c(1, 1, 1, 1), "line")
        )

    g <- ggplot(y.df, aes(x = Hours, y = TPM, color = Gene))
    g <- g + geom_line(aes(group = Gene),size=1.0)
    g <- g + geom_errorbar(aes(ymin = TPM - sd, ymax = TPM + sd, width = 0.3))
    g <- g + theme
    g <- g + guides(color = guide_legend("Gene"), fill = guide_legend("Gene"))
    g <- g + scale_color_jama()
    plot(g)

    ggsave(g, file = paste0(goi, "_.expression.pdf"), width = 9, height = 7)
}