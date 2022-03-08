#===================================================================================================#

# Expression Table

    #For  "Heatmap for methylation regulating TF (Fig. 2b, d)" in Methyl850_HEP.r

#===================================================================================================#
expMatrix <- read.table("iPS_HEP_diff_F5.TPM.txt", header=TRUE, row.name=NULL, sep="\t", stringsAsFactors=FALSE)
expMatrix <- expMatrix[,c(1:4,14:16,5:13)]    # reorder based on time point

#summrize by gene
promoterID <- expMatrix[,"prmtrID"]    # extraction of promoter IDs
gene_name <- sapply(strsplit(promoterID, "@"), function(x){x[2]})    # Extraction of gene names
unique.gene_name <- unique(gene_name)    #removing the redundancy
## Combine all promoter tag counts of each gene
## !!This process takes times!!
gene_expMatrix <- NULL
counter <- 0    # Initialization of the counter
for(i in unique.gene_name){
gene_expMatrix <- rbind(gene_expMatrix, apply(expMatrix[which(gene_name == i),-1],2,sum))
## For loop counter
counter <- counter + 1
progress <- floor((counter/length(unique.gene_name))*100)
cat(paste0(" ##### Progress:",progress, "% #####\r"))
}
rownames(gene_expMatrix) <- unique.gene_name
save(gene_expMatrix, file="gene_expMatrix.RData")    ## aveing the data

#mean and sd
rep_group <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
## Mean
mean.gene_expMatrix <- cbind(apply(gene_expMatrix[,which(rep_group==1)], 1, mean),
    apply(gene_expMatrix[,which(rep_group==2)], 1, mean),
    apply(gene_expMatrix[,which(rep_group==3)], 1, mean),
    apply(gene_expMatrix[,which(rep_group==4)], 1, mean),
    apply(gene_expMatrix[,which(rep_group==5)], 1, mean)
)
colnames(mean.gene_expMatrix) <- c("iPS_HEP_D00", "iPS_HEP_D07", "iPS_HEP_D14", "iPS_HEP_D21", "iPS_HEP_D28")
## SD
sd.gene_expMatrix <- cbind(apply(gene_expMatrix[,which(rep_group==1)], 1, sd),
    apply(gene_expMatrix[,which(rep_group==2)], 1, sd),
    apply(gene_expMatrix[,which(rep_group==3)], 1, sd),
    apply(gene_expMatrix[,which(rep_group==4)], 1, sd),
    apply(gene_expMatrix[,which(rep_group==5)], 1, sd)
)
colnames(sd.gene_expMatrix) <- c("iPS_HEP_D00", "iPS_HEP_D07", "iPS_HEP_D14", "iPS_HEP_D21", "iPS_HEP_D28")

save(list=ls(), file="out/CAGE_HEP/gene_expMatrix.RData")
#===================================================================================================#

# Marker gene exprssion (Supplementary Figure 1c, 2c and 4c)

    #description: This script vidualize expression of indicated genes as a line plot.

#===================================================================================================#
#loding the CAGE expression data
expMatrix <- read.table("iPS_HEP_diff_F5.TPM.txt", header = TRUE, row.name = NULL, sep = "\t", stringsAsFactors = FALSE)
expMatrix <- expMatrix[,c(1:4,14:16,5:13)]    # reorder based on time point

#Compute mean and sd of replicates
mean.expMatrix <- cbind(
    apply(expMatrix[,c(2:4)], 1, mean),
    apply(expMatrix[,c(5:7)], 1, mean),
    apply(expMatrix[,c(8:10)], 1, mean),
    apply(expMatrix[,c(11:13)], 1, mean),
    apply(expMatrix[, c(14:16)], 1, mean)
)
colnames(mean.expMatrix) <- c(0, 7, 14, 21, 28)
rownames(mean.expMatrix) <- expMatrix$prmtrID

#compute sd of replicates
sd.expMatrix <- cbind(
    apply(expMatrix[,c(2:4)], 1, sd),
    apply(expMatrix[,c(5:7)], 1, sd),
    apply(expMatrix[,c(8:10)], 1, sd),
    apply(expMatrix[,c(11:13)], 1, sd),
    apply(expMatrix[, c(14:16)], 1, sd)
)
colnames(sd.expMatrix) <- c(0, 7, 14, 21, 28)
rownames(sd.expMatrix) <- expMatrix$prmtrID

#extraction of GOI
psc_marker <- c("POU5F1", "NANOG")
de_marker <- c("SOX17", "FOXA2")
hep_marker <- c("HNF1A", "TBX3", "HNF4A", "ASGR1", "APOB","SLC10A1", "AFP", "KRT18", "SERPINA1" "ALB", "CYP3A4", "FOXA2")
methylation_genes <- c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3")

gois <- c(psc_marker, de_marker,hep_marker, methylation_genes)

for(i in gois){
    #extraction of GOI expression data
    goi <- i
    gene_names <- sapply(strsplit(rownames(mean.expMatrix), "@"), function(x){ x[2] })

    goi.mean.expMatrix <- mean.expMatrix[gene_names %in% goi, , drop=FALSE]
    goi.sd.expMatrix <- sd.expMatrix[gene_names %in% goi, ,drop=FALSE]

    ##order based on promoter number
    goi.mean.expMatrix <- goi.mean.expMatrix[order(as.numeric(gsub("p", "", sapply(strsplit(rownames(goi.mean.expMatrix), "@"), function(x){x[1]})))),]
    goi.sd.expMatrix <- goi.sd.expMatrix[order(as.numeric(gsub("p", "", sapply(strsplit(rownames(goi.sd.expMatrix), "@"), function(x){x[1]})))),]

    #remove low expreeesion promoters
    cutoff <- 5
    tpm_cutoff <- which((apply(goi.mean.expMatrix, 1, max) >= cutoff))
    goi.mean.expMatrix <- goi.mean.expMatrix[tpm_cutoff, ,drop=FALSE] 
    goi.sd.expMatrix <- goi.sd.expMatrix[tpm_cutoff,, drop = FALSE]

    if(nrow(goi.mean.expMatrix) == 0) next

    #line plot
    goi.mean.expMatrix %>%
        as.data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        pivot_longer(cols = -Gene, values_to = "TPM", names_to = "Days") -> y

    goi.sd.expMatrix %>%
        as.data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        pivot_longer(cols = -Gene, values_to = "sd", names_to = "Days") -> sd.y

    y.df <- cbind(y, sd = sd.y[,3])
    y.df$Days <- factor(y.df$Days, levels = c(0, 7, 14, 21, 28))

    theme <- theme(
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks=element_line(colour="black"),
        legend.key = element_rect(fill = "white"),
        legend.position="bottom",
        plot.margin = unit(c(1, 1, 1, 1), "line"))
        
    g <- ggplot(y.df, aes(x = Days, y = TPM, color = Gene))
    g <- g + geom_line(aes(group = Gene))
    g <- g + geom_errorbar(aes(ymin = TPM - sd, ymax = TPM + sd, width = 0.3))
    g <- g + theme 
    g <- g + guides(color=guide_legend("Gene"),fill=guide_legend("Gene"))
    plot(g)

    ggsave(g, file = paste0(goi, "_.expression.pdf"), width = 9, height = 7)
}