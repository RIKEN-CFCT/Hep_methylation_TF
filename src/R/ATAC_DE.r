#===================================================================================================

# Accesibility change at Demethyl/hyper-methyl maintaining regions (Fig. S6C)

#    peak based analysis
#    1. peak calling: bam2peaks_ATAC.sh (wrapper of MACS2)

#===================================================================================================
#-------------------------------------------------
#overlap peaks
#-------------------------------------------------
peak_dir <- paste0(getwd(), gsub("\\.", "", data_dir), "peakcalling")
bam_dir <- paste0(getwd(), gsub("\\.", "", data_dir), "bam/processed")
gf <- peakreference2(dir = peak_dir, pattern = "filt.narrowPeak")

# Experiment design
experiment_BAMfile <- data.frame(sampleid = c("DE00", "DE48", "DE54", "DE60", "DE66", "DE72"), 
    timepoint =c("0h", "48h", "54h", "60h", "66h", "72h"),
    group = 1:6,
    BAMfile = c("0h.processed.bam", 
        "48h.processed.bam", 
        "54h.processed.bam", 
        "60h.processed.bam", 
        "66h.processed.bam", 
        "72h.processed.bam"))

#-------------------------------------------------
# create a TCA object
#-------------------------------------------------
tca <- TCA(design = experiment_BAMfile, genomicFeature = gf)
tca <- countReads(tca, dir = bam_dir)

#-------------------------------------------------
# convert the TCA object to a SummarizedExperiment object
#-------------------------------------------------
counts <- tca@counts
rowRanges <- GRanges(tca@genomicFeature$chr,
    IRanges(start = tca@genomicFeature$start, tca@genomicFeature$end),
    strand = "*",
    feature_id = tca@genomicFeature$id)
colData <- DataFrame(tca@design)

peak.counts <- SummarizedExperiment(assays = list(counts = counts),
    rowRanges = rowRanges, 
    colData = colData)

# add library size
dgel <- asDGEList(peak.counts)
peak.counts$totals <- dgel$samples$lib.size

# add cpm and log com
cpm <- cpm(dgel, log = FALSE)
assays(peak.counts)$cpm <- cpm

logcpm <- cpm(dgel, log = TRUE)
assays(peak.counts)$logcpm <- logcpm

#-------------------------------------------------
#filtering the peaks opend at 48h
#-------------------------------------------------
opened_peak.counts <- peak.counts[logcpm[,2] - logcpm[,1] >= 2, ]

#-------------------------------------------------
# categolize the opened peaks into 
# demethylated, undemethyl.hyper, and undemethyl.hypo
#-------------------------------------------------
#probe regions
all_probe <- bed2granges("./data/Infinium_anntations/MethylationEPIC_v-1-0_B4.bed")

#demethylated regions at any timepoints 
all_demet <- bed2granges("./data/Methyl850_DE/bed_d1/Demet_all_integrated.bed")

#undemethylated regions
non_demet <- all_probe[!(all_probe$name %in% all_demet$name)]

#undemethylated|ypermethylated regions
selDataMatrix <- read.table("sel_processed_Mval.txt")    #M-value matrix
non_demet_m <- selDataMatrix[rownames(selDataMatrix) %in% non_demet$name, ]
#setting of cut-ff for hyper- nad hypo-methylated probes
non_demet_mean_m <- apply(non_demet_m, 1, mean)
hist(non_demet_mean_m)        #from the histgram, define M > 1 and < - 1 as hyper and hypo methylated regions, respectively

non_demet_hyper <- non_demet[non_demet$name %in% names(non_demet_mean_m[non_demet_mean_m > 1])]
non_demet_hypo <- non_demet[non_demet$name %in% names(non_demet_mean_m[non_demet_mean_m < -2])]

demet_opened_peaks <- opened_peak.counts[findOverlaps(opened_peak.counts@rowRanges, all_demet, ignore.strand = TRUE)  %>% queryHits  %>% unique]
non_demet_hyper_opened_peaks <- opened_peak.counts[findOverlaps(opened_peak.counts@rowRanges, non_demet_hyper, ignore.strand = TRUE)  %>% queryHits %>% unique]

#-------------------------------------------------
#violin plot
#-------------------------------------------------
demet_logcpm <- assays(demet_opened_peaks)$logcpm
colnames(demet_logcpm) <- c("0 h", "48 h", "54 h", "60 h", "66 h", "72 h" )

demet_logcpm  %>% 
    as.data.frame()  %>% 
    mutate(peak_ID = rownames(demet_logcpm))  %>% 
    pivot_longer(-peak_ID, names_to = "timepoint", values_to = "log_CPM") -> demet_logcpm_df

p1 <- ggplot(demet_logcpm_df, aes(x = timepoint, y = log_CPM, colour= timepoint))
p1 <- p1 + geom_violin(trim = T, fill = "#999999", linetype = "blank", alpha = I(1/3))
p1 <- p1 + stat_summary(geom = "pointrange", 
    fun = mean, 
    fun.min = function(x){
        mean(x)-sd(x)},
    fun.max = function(x){
        mean(x)+sd(x)},
    size=1,
    alpha=.5)
p1 <- p1 + stat_summary(fun = mean, geom="line", aes(group=1))
p1 <- p1 + ylim(-2.5,5.5)
p1 <- p1 + my.theme + theme(legend.position = 'none')
p1 <- p1 + 	scale_color_jama()
p1 <- p1 + labs(title="Demethylated regions", y="log CPM")

hyper_non_demet_logcpm<- assays(non_demet_hyper_opened_peaks)$logcpm
colnames(hyper_non_demet_logcpm) <- c("0 h", "48 h", "54 h", "60 h", "66 h", "72 h" )

hyper_non_demet_logcpm  %>% 
    as.data.frame()  %>% 
    mutate(peak_ID = rownames(hyper_non_demet_logcpm))  %>% 
    pivot_longer(-peak_ID, names_to = "timepoint", values_to = "log_CPM") -> hyper_non_demet_logcpm_df

p2 <- ggplot(hyper_non_demet_logcpm_df, aes(x = timepoint, y = log_CPM, colour= timepoint))
p2 <- p2 + geom_violin(trim = T, fill = "#999999", linetype = "blank", alpha = I(1/3))
p2 <- p2 + stat_summary(geom = "pointrange", 
    fun = mean, 
    fun.min = function(x){
        mean(x)-sd(x)},
    fun.max = function(x){
        mean(x)+sd(x)},
    size=1,
    alpha=.5)
p2 <- p2 + stat_summary(fun = mean, geom="line", aes(group=1))
p2 <- p2 + ylim(-2.4,5.5)
p2 <- p2 + my.theme + theme(legend.position = 'none', axis.title.y = element_blank())
p2 <- p2 + 	scale_color_jama()
p2 <- p2 + labs(title="Hyper-methylated regions")

p <- plot_grid(p1, p2, nrow = 1, align = "hv",axis = "tblr")
print(p)
