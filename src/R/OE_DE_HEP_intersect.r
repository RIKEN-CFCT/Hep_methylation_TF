
#===================================================================================================

# intersect (Supplementary Figure 6e)

#===================================================================================================
#methylation array data loading
oe_data <- read.table("Methyl850_OE_sel_processed_Mval.txt")    #M-value matrix of GATA6 overexpression 
hep7_data <- read.table("Methyl850_HEP_sel_processed_Mval.txt")    #M-value matrix of hepatocyte differentiation time-course

#common probes
oe_data <- oe_data[rownames(oe_data) %in% rownames(hep7_data),]
hep7_data <- hep7_data[rownames(hep7_data) %in% rownames(oe_data),]

#all probes
all_probe_posi <- na.omit(probeID2position(rownames(oe_data)))
all_probe_gr <- GRanges(seqnames = all_probe_posi[, 1],
                        ranges = IRanges(start = all_probe_posi[, 2],
                                         end = all_probe_posi[, 2] + 1),
                        strand = "*",
                        name = rownames(all_probe_posi))
all_probe_gr <- all_probe_gr + 200

#demethylated probes in OE
oe_demet_id <- DmpId(selDataMatrix = oe_data, ControlColnum = 1, TreatmentColnum = 2)
oe_posi <- probeID2position(oe_demet_id)
oe_gr <- GRanges(seqnames = oe_posi[,1],
                        ranges = IRanges(start = oe_posi[,2],
                                         end = oe_posi[,2] + 1),
                        strand = "*",
                        name = rownames(oe_posi))
oe_gr <- oe_gr + 200

#demethylated probes at day 7 comp to day 0
hep7_demet_id <- DmpId(selDataMatrix = hep7_data, ControlColnum = 1, TreatmentColnum = 2)
hep7_posi <- probeID2position(hep7_demet_id)
hep7_gr <- GRanges(seqnames = hep7_posi[,1],
                        ranges = IRanges(start = hep7_posi[,2],
                                         end = hep7_posi[,2] + 1),
                        strand = "*",
                        name = rownames(hep7_posi))
hep7_gr <- hep7_gr + 200

#overlaped demethylated probes
intersect_posi <- probeID2position(oe_demet_id[oe_demet_id %in% hep7_demet_id])
intersect_gr <- GRanges(seqnames = intersect_posi[,1],
                        ranges = IRanges(start = intersect_posi[,2],
                                         end = intersect_posi[,2] + 1),
                        strand = "*",
                        name = rownames(intersect_posi))
intersect_gr <- intersect_gr + 200

#GATA6 ChIP
gata6_gr <- bed2granges("GATA6_72h_peaks.narrowPeak")

#-------------------------------------------------
#venn's diagrams
#-------------------------------------------------
#compute overlaps
(n.oe <- length(oe_gr))
(n.hep7 <- length(hep7_gr))
(n.GATA6 <- length(gata6_gr))
(nOE_hep7 <- length(intersect_gr))
(nhep7_GATA6 <- length(findOverlaps(gata6_gr, hep7_gr)))
(nGATA6_OE <- length(findOverlaps(gata6_gr, oe_gr)))
(n.all <- length(findOverlaps(gata6_gr, intersect_gr)))

#OE vs GATA6
v <- venneuler(c(OE = n.oe, GATA6 = n.GATA6, "GATA6&OE" = nGATA6_OE))
pdf(paste0(out_dir, "OE_GATA6_intersect.pdf"))
plot(v, main = "Overlaps of demethylated regions", col = c("deeppink4", "azure4"))
dev.off()

#DE vs GATA6
v <- venneuler(c(DE = n.hep7, GATA6 = n.GATA6, "DE&GATA6" = nhep7_GATA6))
pdf(paste0(out_dir, "DE_GATA6_intersect.pdf"))
plot(v, main = "Overlaps of demethylated regions", col = c("blue", "green"))
dev.off()

#DE_GATA6 vs OE
v <- venneuler(c(DE_GATA6 = nhep7_GATA6, OE = n.oe, "DE_GATA6&OE" = n.all))
pdf(paste0(out_dir, "DEGATA6_OE_intersect.pdf"))
plot(v, main = "Overlaps of demethylated regions", col = c("seagreen", "royalblue3"))
dev.off()

#Three intersects
v <- venneuler(c(OE = n.oe, DE = n.hep7, GATA6 = n.GATA6, "OE&DE" = nOE_hep7, "DE&GATA6" = nhep7_GATA6, "GATA6&OE" = nGATA6_OE, "OE&DE&GATA6" = n.all))
pdf(paste0(out_dir, "OE_DE_GATA6_intersect.pdf"))
plot(v, main = "Overlaps of demethylated regions", col = c("blue", "green", "orange"))
dev.off()

#-------------------------------------------------
#bed files
#-------------------------------------------------
all_gr <- intersect_gr[subjectHits(findOverlaps(gata6_gr, intersect_gr))]
export.bed(all_gr, con=paste0(data_dir, 'all_intersects.bed'))

#annotation
anno <- read.csv("./data/MethylationEPIC_v-1-0_B4.csv", skip = 7, header=TRUE)
all_anno <- anno[anno[,1] %in% all_gr$name, c(1,11,12,13,15, 16, 17,18, 19,20,21, 22)]
write.table(all_anno, file = paste0(out_dir, "all_intersect_anno.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

#-------------------------------------------------
#permutation test
#-------------------------------------------------
pt1 <- permTest(A = hep7_gr, 
    B = gata6_gr, 
    universe = all_probe_gr,
    randomize.function = resampleRegions,
    evaluate.function = numOverlaps,
    ntimes = 1000,
    verbose = TRUE)
pt1

GATA6_DE_intersect_gr <-hep7_gr[subjectHits(findOverlaps(gata6_gr, hep7_gr))]
pt2 <- permTest(A = GATA6_DE_intersect_gr, 
    B = oe_gr, 
    universe = all_probe_gr,
    randomize.function = resampleRegions,
    evaluate.function = numOverlaps,
    ntimes = 1000,
    verbose = TRUE)
pt2
