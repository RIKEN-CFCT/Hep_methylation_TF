#-----------------------------------------
#DMPrGREAT
#-----------------------------------------
#' This function peform ontology analysis based on GREAT and retreaves GOs and panther pathway
#' 
#' This function peform ontology analysis based on GREAT for differentially methylated probes of infinium methylation arrays. The script outputs enriched GOs and panther pathway table text files to a indicated and returns resulting rGREAT object.
#' 
#' @param data a M-value data matrix
#' @param out_prefix output directory putprefix
#' @param seq_range ditance from Cytocine to make bed file
#' @param ControlColnum the control column of M-value matrix 
#' @param TreatmentColnum the treatment column of M-value matrix
#' @param version version of Infinium Methjylation array (450 or 850(orEPIC))
#' @param sampling If sampling number is indicated and DMP is more than sampling number, analysis is ran with randomly selected DMP. If FALSE, all DMPs are analyzed.
#' 
#' @importFrom InfiniumDiffMetMotR DmpId probeID2position
#' @importFrom dplyr %>% select
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rGREAT submitGreatJob getEnrichmentTables getEnrichmentTables
#' 
#' @return rGREAT result
#' @keywords GO, GREAT, DMP
#' @export
#' 

DMPrGREAT <- function(data, out_prefix="GREAT_result", ControlColnum, TreatmentColnum, version = "850", sampling = FALSE, seq_range = c(-100, 100),...){
    #DMP identification
    cat("DMP identification...\n")
    selDataMatrix <- data
    DMP_IDs <- DmpId(selDataMatrix = selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum,...)

    #sammpling
    if(!sampling == FALSE & length(DMP_IDs) >= sampling){ 
        DMP_IDs <- DMP_IDs[floor(runif(sampling,1,length(DMP_IDs)))]
    }

    #convert ID to position
    allProbe_IDs <- rownames(selDataMatrix)
    if(version=="450"){
        target_position <- na.omit(probeID2position(probe_IDs=DMP_IDs,anno_info=Methyl450anno))	#conversion of DMP IDs to position
        anno_info <- Methyl450anno
        CHR37 <- unlist(anno_info %>% select("CHR"))
        CHR37 <- paste("chr", CHR37, sep="")
        CPG37 <- anno_info %>% select("MAPINFO")
        bg_position <- na.omit(cbind(CHR37, CPG37))
        positionsList <- list("target" = target_position, "background" = bg_position)    #integrate DMP positoins and NC positions
    }else if ((version=="EPIC")||(version=="850")){
        target_position <- na.omit(probeID2position(probe_IDs=DMP_IDs, anno_info=EPICanno))	#conversion of DMP IDs to position
        anno_info <- EPICanno
        CHR37 <- unlist(anno_info %>% select("CHR"))
        CHR37 <- paste("chr", CHR37, sep="")
        CPG37 <- anno_info %>% select("MAPINFO")
        bg_position <- na.omit(cbind(CHR37, CPG37))
        positionsList <- list("target" = target_position,  "background" = bg_position)	#integrate DMP positoins and NC positions
    }

    bed_list <- lapply(positionsList, function(x){position2bed3(positions = x, seq_range = seq_range)})

    #convert position to GRanges
    target_bed_gr <- with(bed_list[[1]], GRanges(chr, IRanges(start,end)))
    bg_bed_gr <- with(bed_list[[2]], GRanges(chr, IRanges(start,end)))

    ##GREAT analysis
    cat("Running GREAT...\n")
    job <- submitGreatJob(target_bed_gr, bg = bg_bed_gr)
    go_tb = getEnrichmentTables(job)
    mp_tb = getEnrichmentTables(job, ontology = "Mouse Phenotype")

    ## output of data
    cat("Result writing...\n")
    dir.create(out_prefix, recursive = TRUE)

    go_mf_out <- paste0(out_prefix, "/GO_MF_table.txt")
    go_bp_out <- paste0(out_prefix, "/GO_BP_table.txt")
    go_cc_out <- paste0(out_prefix, "/GO_CC_table.txt")
    mp_out <- paste0(out_prefix, "/Mouse_Phenotype_table.txt")

    write.table(go_tb[["GO Molecular Function"]], file=go_mf_out, sep="\t", quote=FALSE, row.names = FALSE)
    write.table(go_tb[["GO Biological Process"]], file=go_bp_out, sep="\t", quote=FALSE, row.names = FALSE)
    write.table(go_tb[["GO Cellular Component"]], file=go_cc_out, sep="\t", quote=FALSE, row.names = FALSE)
    write.table(mp_tb, file=mp_out, sep="\t", quote=FALSE, row.names = FALSE)

    cat("Analysis Succeeded!\n")
    return(job)
}

#--------------------------------------------
# dmr2bed
#--------------------------------------------
#' Reading cellranger processed 10x data as a seurat object.
#' 
#' Read the cell ranger processed 10x data, assign a sample name, normalize, and find veriable genes
#' 
#' @param data_amatrix path of ceranger preocessed result (one upper from "outs" directory)
#' @param cutoff cutoff for delta-M-value
#' @param p.cutoff cutoff for Welch's t-test
#' @param ControlColnum column number for Control sample(s)
#' @param TreatmentColnum  column number for target sample(s)
#' @param MethylDemethyl direction of Ddifferentially methylated regions (methylated or demethylated)
#' @param version version of Infinium Methjylation array (450 or 850(orEPIC))
#' @param extension a vector containing dounstream and upsteream extensions
#' 
#' @importFrom InfiniumDiffMetMotR DmpId probeID2position
#' 
#' @return a bed format-like data.frame
#' 
#' @keywords DMR bed
#' @export 

dmr2bed	<- function(data_matrix, cutoff = 2, p.cutoff = 0.001, ControlColnum, TreatmentColnum, MethylDemethyl="Demethyl", version = "850", extension = c(-0, 1)){
    ## Identification of DMR
	DMP_IDs <- DmpId(selDataMatrix=data_matrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, p.cutoff=p.cutoff, cutoff= cutoff, MethylDemethyl=MethylDemethyl)

    ## Convert ID to position
	if(version=="450"){
		target_position <- na.omit(probeID2position(probe_IDs=DMP_IDs,anno_info=Methyl450anno))	#conversion of DMP IDs to position
	}else if ((version=="EPIC")||(version=="850")){
		target_position <- na.omit(probeID2position(probe_IDs=DMP_IDs, anno_info=EPICanno))	#conversion of DMP IDs to position
	}
    
    ## dM computation
    data_matrix.dmr <-  data_matrix[rownames(data_matrix) %in% DMP_IDs,]
    dm <- data_matrix.dmr[,TreatmentColnum] - data_matrix.dmr[,ControlColnum]
    
    ## bed formatting
	CHR37 <- target_position[,1]
	Probe_posi <- as.numeric(target_position[,2])
	Str37 <- Probe_posi + extension[1]
	End37 <- Probe_posi + extension[2]
	M.val <- dm
	strand <- rep("*", length(dm))

	bed <- data.frame(CHR37, Str37, End37, DMP_IDs, M.val, strand)
	out <- bed
	return(out)
}

#--------------------------------------------
# tcHeatmap
#--------------------------------------------
#' timecourse value data frame heatmap
#' 
#' This function generate a heatmap based on the timecourse value data frame 
#' 
#' @param df a data frame of m-calue
#' @param mycols color pallet for the heatmap. defaut is PuBu
#' @param z burbose TRUE convert the value to z-score
#' 
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid unit gpar
#' 
#' @return a complexHeatmap object
#' 
#' @keywords M-value heatmap
#' @export
#' 
#

tcHeatmap <- function(df, mycols, z = TRUE,...){
	if(z == TRUE){
		df <- genescale(df, axis = 1, method = "Z")
		name <- "M (z-score)"
	}else{
		df <- df
		name <- "M-value"
	}

    hc.rows <- hclust(dist(as.matrix(replace((df), is.na(df), 0)))) # クラスタリング

    ht <- Heatmap(df,
        width = unit(10, "cm"),    #With of heatmap
        name = name,    #title of legend
        column_title = "Timepoint",
        column_title_gp = gpar(fontsize = 15),
        row_title = "Probes",
        row_title_gp = gpar(fontsize = 15),
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 15),    #Text size for row names
        column_names_gp = gpar(fontsize = 15),    #Text size for column names
        col = mycols,
        cluster_rows = hc.rows,
        show_row_dend = FALSE,
        row_dend_reorder = FALSE,
        cluster_columns = FALSE,
        heatmap_legend_param = list(
            title = name,
            title_gp = gpar(fontsize = 15),
            title_position = "topleft",
            labels_gp = gpar(fontsize = 15),
            legend_height = unit(5, "cm"),
            grid_width = unit(1, "cm")
        ),
		...
    )

	return(ht)
}


#--------------------------------------------
# tcLinePlot
#--------------------------------------------
#' timecourse averaged line plot
#' 
#' This function generate a averaged lienplot based on the timecourse value data frame 
#' 
#' @param df a data frame of m-calue
#' @param color color of line
#' @param fill color of sd shade
#' 
#' @importFrom ggplot2
#' 	ggplot
#' 	geom_ribbon
#' 	aes
#' 	geome_line
#' 	geome_point
#' 	labs
#' 	scale_color_manual
#' 	scale_fill_manual
#' 	ylim
#' 	theme
#' 	element_text
#' 	element_line
#' 	element_rect
#' 	ggsave
#' @importFrom dplyr %>% mutate
#' 
#' @return a ggplot object
#' 
#' @keywords M-value line plot
#' @export
#' 
#'     
tcLinePlot <- function(df, is.log = FALSE, color, fill, y_lim = FALSE){
    df1 <- df + 0.0001

	#convert to fold-change (vs most left column)
		df2 <- NULL
	if(is.log == TRUE){
		for (j in seq(nrow(df1))){
    	    norm <- df1[j,]-df1[j,1]
        	df2 <- rbind(df2, norm)
    	}
	}else{
		for (j in seq(nrow(df1))){
			norm <- log2(df1[j,]/df1[j,1])
			df2 <- rbind(df2, norm)
		}
	}
    rownames(df2) <- rownames(df1)
    colnames(df2) <- colnames(df1)

	#compute mean and sd
    df3 <- data.frame(
        timepoint = as.numeric(colnames(df2)),
        mean = sapply(df2, mean), 
        sd = sapply(df2, sd) 
    )
    df3 <- mutate(df3, ymin=df3$mean-df3$sd) %>% 
        mutate(ymax=mean+sd)

    g <- ggplot(df3, aes(x = timepoint))
    g <- g + geom_ribbon(aes(ymin = ymin, ymax=ymax), fill= fill, color=NA)
    g <- g + geom_line(aes(y=mean), color=color, size=1)
    g <- g + geom_point(aes(y=mean), color="black", size=1.7)
    g <- g + labs(x = "Timepoint", y = "Fold-Change (log)")
	if(y_lim != FALSE){
    	g <- g + ylim(y_lim[1], y_lim[2])
	}
    g <- g + theme(
        axis.text = element_text(color = "black", size = 18),    # axis texit setting
        axis.ticks = element_line(color = "black", size = .5),    # Setting of scale ticks line
        axis.title = element_text(color = "black", size = 18),    # Setting of scale ticks texit
        panel.background = element_rect(fill = "white", color = NA),    # panel background color
        panel.border = element_rect(fill = NA, color = "black"),    # panel border
        panel.grid.major = element_line(color = "gray90"),    # panel major grid
        panel.grid.minor = element_line(color = "gray90"),   # panel minor grid 
        legend.position = "none"
    )

	return(g)
}

#--------------------------------------------
# covHeatmap
#--------------------------------------------
#' heatmap of coverage distribution around given target reagions
#' 
#' This function make a complexHeatmap object that is a heatmap of coverage (read) distribution. This funciton is awrapper of enrichedHeatmap.
#' 
#' @param target a path of target bed file
#' @param bw a path of bigWig coverage file
#' @param extension a numeric. up- and down-stream ranges of heatmap
#' @param col_fun color function for the heatmap
#' @param add_profile a logic wheather add frofile plot at upper of the heatmap
#' 
#' @importFrom data.table fread
#' @importFrom rtracklayer BigWigSelection import
#' @importFrom GenomicRanges GRanges resize
#' @importFrom IRanges IRanges
#' @importFrom EnrichedHeatmap normalizeToMatrix EnrichedHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom complexHeatmap HeatmapAnnotation draw
#' 
#' @return a GRanges object
#' 
#' @keywords bed GRanges
#' @export 
#' 
#' 

covHeatmap <- function(
    target,
    bw,
    extension = 5000,
    cols = c("darkblue", "darkgoldenrod1"),
    add_profile = TRUE, ...
    ){
    library("EnrichedHeatmap")
    library("data.table")
    library("rtracklayer")
    library("circlize")

    #loading of a target file
    tmp.targets <- makeGRangesFromDataFrame(
        df = fread(target,
            stringsAsFactors = F,
            sep = "\t",
            header = F
        ),
        seqnames.field = "V1",
        start.field = "V2",
        end.field = "V3"
    )

    #analysis range
    tmp.extension <-  extension

    #extend center of the targets
    tmp.targets_extended <- resize(tmp.targets, fix = "center", width = tmp.extension * 2)

    #load the bigwig file
    tmp.bigwig <- rtracklayer::import(bw,
                                        format = "BigWig",
                                        selection = BigWigSelection(tmp.targets_extended))

    normMatrix <- normalizeToMatrix(
        signal = tmp.bigwig,
        target = resize(tmp.targets, fix = "center", width = 1),
        background = 0,
        keep = c(0, 0.99), ## minimal value to the 99th percentile
        target_ratio = 0,
        mean_mode = "w0", ## see ?EnrichedHeatmap on other options
        value_column = "score", ## = the name of the 4th column of the bigwig
        extend = tmp.extension
    )

    #heatmap
    #coror of the heatmap
    col_fun = circlize::colorRamp2(quantile(normMatrix, c(0, .99)), cols)

        enrHtmp <- EnrichedHeatmap(
        mat = normMatrix,
        pos_line = FALSE, ## no dashed lines around the start
        border = FALSE, ## no box around heatmap
        col = col_fun, ## color gradients from above

        ## turn off background colors
        rect_gp = gpar(col = "transparent"),

        if (add_profile == TRUE){
            ## options for the profile plot on top
            top_annotation = HeatmapAnnotation(
                enriched = anno_enriched(
                gp = gpar(col = "black", lty = 1, lwd=2),
                col = "black"
                )
            )
        },

        #output setting
        use_raster = TRUE,
        raster_quality = 5,
        raster_device = "png", ...
    )
    return(enrHtmp)
}


#--------------------------------------------
# bedExport
#--------------------------------------------
#' export a bed-like table to a bed file
#' 
#' This function exports a bed-like table to a bed file.
#' 
#' @param object a bed-like table containing chromosome, start, end, ID, and value
#' @param file a name of utput bed file
#' 
#' @return bed file
#' 
#' @keywords bed
#' @export

bedExport <- function(object, file){
    write.table(
        x = object,
        file = file,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
    )
}

#--------------------------------------------
# peakreference2

    #modified peakreference of TCseq. Thanks to Mengjun Wu, Lei Gu
#--------------------------------------------
#' combine and merge multiple BED files
#'
#' This function merges genomic coordinates of a given data frame or
#' reads in BED files (e.g. generated from a peak caller) under given
#' directory and merge genomic regions that have overlapping genomic
#' intervals into a single feature. The single feature represents
#' the widest genomic interval that covers all merged regions.
#'
#' @param data a data frame containg coordinates information of peaks
#' to be merged. Columns of the data frame should be consistent with
#' the BED format where the first column is the name of the chromosome,
#' the second column is the starting position and the third column is
#' the ending position.
#'
#' @param dir character string giving the directory where BED files
#' are stored. If \code{data} is not given, the function will reads
#' in the BED files under \code{code}.
#'
#' @param pattern an \code{\link{regular expression}}, only files that
#' have names match the regular expression will be read in.
#'
#' @param merge logical indicating whether to merge overlapped regions
#' or not. If False, regions are simply combined.
#'
#' @param overlap a numberic value giving the least number of base(s)
#' two regions should overlap when merging them.
#'
#' @param ratio a numberic value giving the thresold of overlapping
#' ratio between two regions to merge them. See '\code{Details}' below
#' for the definition of the overlapping ratio.
#'
#' @return a data frame with four columns: \code{chr}, \code{start},
#' \code{stop}, \code{id}
#'
#' @details
#' The overlapping ratio (OR) is defined as:
#'
#' \deqn{ OR = \frac{n}{\min(length(a), length(b)}}
#'
#' \eqn{a}, \eqn{b} are two genomic regions, \eqn{n} is the number of
#' overlapping bases between region \eqn{a} and region \eqn{b}.
#'
#' @export

peakreference2 <- function (data = NULL, dir = NULL, pattern = NULL, merge = TRUE, 
    overlap = 1, ratio = NULL) 
{
    if (is.null(data) && is.null(dir)) {
        stop("Either a data.frame of all peak information or directory where the BED files store should be given")
    }
    if (!is.null(data)) {
        checkBEDformat(data)
        data[, 1] <- factor(data[, 1])
        peakset <- data
    }
    if (is.null(data) && !is.null(dir)) {
        old <- setwd(tempdir())
        on.exit(setwd(old), add = TRUE)
        setwd(dir)
        filenames <- list.files(pattern = pattern)
        if (length(filenames) == 0) {
            err <- paste0("Can not find file names containing '", 
                pattern, "'.")
            stop(err)
        }
        datalist <- lapply(filenames, function(x) {
            read.table(file = x, header = FALSE)
        })
        peakset <- do.call(rbind, datalist)
        checkBEDformat(peakset)
    }
    peakset <- peakset[order(peakset[, 1], peakset[, 2]), ]
    if (merge) {
        if (overlap <= 0 || round(overlap) != overlap) {
            stop("\"overlap\" must be integer and greater than 0.")
        }
        peakset.sub <- split(peakset, peakset[, 1], drop = TRUE)
        
        #level <- levels(peakset[, 1])
        level <- names(peakset.sub)    #bugfix from peakreference
        
        mergedpeak <- c()
        for (i in seq_len(length(peakset.sub))) {
            temp <- peakset.sub[[i]]
            if (is.null(ratio)) {
                submerge <- intervalmerge(temp[, 2], temp[, 3], 
                  overlap = overlap)
            }
            else {
                submerge <- intervalmerge(temp[, 2], temp[, 3], 
                  ratio = ratio)
            }
            chr <- rep(level[i], length(submerge[, 1]))
            submerge1 <- data.frame(chr, submerge)
            mergedpeak <- rbind(mergedpeak, submerge1)
        }
        name <- paste0("peak", seq_len(length(mergedpeak[, 1])))
        mergedpeak <- data.frame(mergedpeak, name)
        colnames(mergedpeak) <- c("chr", "start", "end", "id")
        mergedpeak
    }
    else {
        peakset
    }
}

checkBEDformat <- function(data){
    if(ncol(data) < 3){
        stop("At least three columns should be provided. The first column contains chromosome name,
         the second column contains starting position, the third column contains ending position.")
    }
    if(class(as.vector(data[, 1])) != "character"){
        stop("The first column contains chromosome name and must be character.")
    }
    if(any(round(data[, 2]) != data[, 2]) &&
        any(round(data[, 3]) != data[, 3])){
        stop("the second and third column contain starting and ending positions, must be numeric.")
    }
}

intervalmerge <- function(a0, b0, overlap = NULL,
                          ratio = NULL){
    if(length(a0) > 1){
        a1 <- c(a0[1])
        b1 <- c(b0[1])
        merge <- NULL
        for(i in seq_len(length(a0) - 1)){
            if(is.null(ratio)){
                if(b1[length(b1)] - a0[i + 1] < overlap){
                    a1 <- append(a1, a0[i + 1])
                    b1 <- append(b1, b0[i + 1])
                } else{
                    b1[length(b1)] <- max(b1[length(b1)], b0[i + 1])
                }
            }
            if(is.null(overlap)){
                len <- min((b1[length(b1)] - a1[length(b1)]),
                   (b0[i + 1] - a0[i + 1]))
                rt <- (b1[length(b1)] - a0[i + 1]) / len
                if(rt < ratio){
                    a1 <- append(a1, a0[i + 1])
                    b1 <- append(b1, b0[i + 1])
                } else{
                    b1[length(b1)] <- max(b1[length(b1)], b0[i + 1])
                }
            }
        }
        merge <- cbind(a1, b1)
    }
    if(length(a0) <= 1){
        a1 <- c(a0[1])
        b1 <- c(b0[1])
        merge <- cbind(a1, b1)
    }
    merge
}

#--------------------------------------------
# kclust
#--------------------------------------------
#' k-mean clustering
#' 
#' This function performs k-mean clustering and plot the results
#' 
#' @param data a matrix or data.frame
#' @param knum a numeric. number of cluster
#' @param out_dir a charactor. directory for output
#' 
#' @importFrom stats kmeans
#' @importFrom tidyverse %>% mutate pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line xlab ylab labs scale_colour_gradientn
#' @importFrom ggsci scale_color_aaas
#' 
#' @return alist of k-mean clustering results, pdfs of expressin patten of each cluster, a pdf showning average expression pattern
#' 
#' @keywords k-means
#' @export 

kclust <- function(data, knum = 3, out_dir = "mean_kcluster/"){
    #create directry
    if(!file.exists(out_dir)){
        dir.create(out_dir, recursive = TRUE)
    }
    
    #clustering
    set.seed(20)
    kClust <- kmeans(data, centers = knum, nstart = 1000, iter.max = 20)
    kClusters <- kClust$cluster

    clust.centroid = function(i, dat, clusters){
        ind = (clusters == i)
        colMeans(dat[ind,])
    }
    kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, data, kClusters)

    #get in long form for plotting
    as.data.frame(kClustcentroids) %>%
        mutate(sample = rownames(kClustcentroids)) %>%
        pivot_longer(co = -sample, names_to = "cluster", values_to = "value") -> Kmolten

    #plot average expression pattern of each k-means cluster
    theme <- theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
    )

    p1 <- ggplot(Kmolten, aes(x = sample, y = value, group = cluster, colour = as.factor(cluster)))
    p1 <- p1 + geom_point()
    p1 <- p1 + geom_line()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("cpm")
    p1 <- p1 + labs(title = "Cluster CPM by Time", color = "Cluster")
    p1 <- p1 + theme + scale_color_aaas()
    p1
    ggsave(p1, file = paste0(out_dir, "kclust_exp_by_time.pdf"))

    kclust_data_list <- NULL
    for(i in seq(knum)){
        #Subset a k-means cluster}
        core <- Kmolten[Kmolten$cluster == i,]

        #get cluster data
        kclust_data <- (data[kClusters == i,])
        #calculate the correlation with the core
        score <- apply(kclust_data, 1, function(x){ cor(x, core$value) })
        #get the data frame into long format for plotting
        as.data.frame(kclust_data) %>%
            mutate(gene = rownames(kclust_data)) %>%
            pivot_longer(co = -gene, names_to = "cluster", values_to = "value") -> kclust_data_molten
        #add the score
        kclust_data_molten <- merge(kclust_data_molten, score, by.x = 'gene', by.y = 'row.names', all.x = T)
        colnames(kclust_data_molten) <- c('gene', 'sample', 'value', 'score')
        #order the dataframe by score
        #to do this first create an ordering factor
        kclust_data_molten$order_factor <- 1:length(kclust_data_molten$gene)
        #order the dataframe by score
        kclust_data_molten <- kclust_data_molten[order(kclust_data_molten$score),]
        #set the order by setting the factors
        kclust_data_molten$order_factor <- factor(kclust_data_molten$order_factor, levels = kclust_data_molten$order_factor)

        # Everything on the same plot
        p2 <- ggplot(kclust_data_molten, aes(x = sample, y = value))
        p2 <- p2 + geom_line(aes(colour = score, group = gene))
        p2 <- p2 + scale_colour_gradientn(colours = c('blue1', 'red2'))
        p2 <- p2 + geom_line(data = core, aes(sample, value, group = cluster), color = "black", inherit.aes = FALSE) #add the core 
        p2 <- p2 + xlab("Time")
        p2 <- p2 + ylab("Expression")
        p2 <- p2 + labs(title = paste0("Cluster ", i , " Expression"), color = "Score")
        p2 <- p2 + theme_bw()
        p2
        ggsave(p2, file = paste0(out_dir, "K", i, "_exp.pdf"))

        kclust_data_list <- c(kclust_data_list, list(kclust_data))
    }
    names(kclust_data_list) <- paste0("K", seq(knum))
    return(kclust_data_list)
}

#--------------------------------------------
# GATA6_motif_analysis
#--------------------------------------------
#' Analyze GATA6 motif overrepersentation and plot
#' 
#' This script compute enrichment score of GATA6 for the DMR and returen ggplot object of enrichment score
#' 
#' @param selDataMatrix data Matrix
#' @param prefix prefix of output files
#' @param MD "Methyl" or "Demethyl". indication of methylated DMR or demethylated DMR
#' 
#' @importFrom InfiniumDiffMetMotR DmpId motifDBList probeID2position stratSampling seqExtract writeSplitSeq splitSeqMotDist enrichScoreDist
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr  %<>% 
#' @importFrom ggplot2 ggplot aes geom_line ylim xlab ylab geom_hline theme
#' 
#' @return a ggplot object
#' 
#' @keywords GATA6 otif overrepresentation
#' @export 

GATA6_motif_analysis <- function(selDataMatrix, prefix, MD = c("Methyl", "Demethyl")){
    outname <- prefix
    #-------------------------------------------------
    #motif database construction
    #-------------------------------------------------
    motif_names <- "GATA6"
    motifDBList <- IMAGE_PWMlist

    motifDBList <- motifDBList[grep(motif_names, names(motifDBList))]

    #-------------------------------------------------
    #extraction of differentially methylated probe ID
    #-------------------------------------------------
    DMP_IDs <- DmpId(selDataMatrix, ControlColnum = 1, TreatmentColnum = 2, MethylDemethyl = MD)
    nDMP_IDs <- length(DMP_IDs)
    allProbe_IDs <- rownames(selDataMatrix)

    #-------------------------------------------------
    #extraction of DMP positions and stratified sampling
    #-------------------------------------------------
    probe_annotation <- InfiniumDiffMetMotR::EPICanno
    target_position <- na.omit(probeID2position(probe_IDs = DMP_IDs, anno_info = probe_annotation)) 
    randomProbe_IDs <- stratSampling(target_IDs = DMP_IDs, anno_info = probe_annotation) 
    random_position <- na.omit(probeID2position(probe_IDs = randomProbe_IDs, anno_info = probe_annotation)) 
    positionsList <- list("target" = target_position, "random" = random_position)
    ##write DMP position
    DMPtOutfile <- paste0(out_dir, outname, '_DMP_position.txt')
    out_target_position <- cbind(probeID = rownames(target_position), target_position)
    write.table(out_target_position, file = DMPtOutfile, quote = FALSE, row.names = FALSE, col.names = TRUE)

    #-------------------------------------------------
    #Sequence extractio
    #-------------------------------------------------
    #Read human hg19 genomic sequence
    tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep = ":"))
    genome <- eval(parse(text = tmp))

    #sequence extraction
    seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
    sequences <- lapply(positionsList, function(x){ seqExtract(positions = x, genome = genome, seq_range) })

    #-------------------------------------------------
    #Write splitted sequence
    #-------------------------------------------------
    #make a templrary outpput directory
    tempDir <- paste(outname, "_temp", sep = "")
    dir.create(tempDir)

    #writing the sequences to splitted files
    seqs <- sequences$target
    target_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "target")
    seqs <- sequences$random
    random_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "random")

    rm(sequences)
    rm(seqs)
    invisible(replicate(3, gc()))
    #-------------------------------------------------
    #motif search
    #-------------------------------------------------
    target_positionsList <- splitSeqMotDist(filenames = target_all_filenames, motif_list = motifDBList, min.score = "90%")
    file.remove(target_all_filenames)
    gc()
    random_positionsList <- splitSeqMotDist(filenames = random_all_filenames, motif_list = motifDBList, min.score = "90%")
    file.remove(random_all_filenames)
    gc()

    #-------------------------------------------------
    #computation of enrichment score
    #-------------------------------------------------
    #calculation
    motif_name <- names(motifDBList)
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

    #plot
    ranks <- seq(seq_range[1], seq_range[2], length=1001)
    data <- as.data.frame(enrichment_scores)
    colnames(data) <- "ES"
    data <- mutate(data, dist = ranks)

    data  %<>% 
        pivot_longer(-dist, names_to = "sample", values_to = "ES")

    g <- ggplot(data, aes(x = dist, y = ES,  group = sample, colour = sample))
    g <- g + geom_line()
    g <- g + ylim(c(-0.04,0.5))
    g <- g + xlab("Distance from CpG")
    g <- g + ylab("Enrichment Score")
    g <- g + geom_hline(yintercept=0)
    g <- g + my.theme + theme(legend.position= 'none')

    file.remove(tempDir)
    return(g)
}
