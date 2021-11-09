#===================================================================================================#

#Preprocessing

#===================================================================================================#
selDataMatrix <- lumiMethyNorm(
    fileName = "Matrix_signal_intensities.txt",
    inputtype = "signal",
    sample_names = c("mock", "GATA6")
)
#===================================================================================================#

#Scatter plot (Fig. 3B)

#===================================================================================================#



treatment_data <-  selDataMatrix[, 2]
control_data <- selDataMatrix[, 1]
main <- "GATA6 Overexpression"
xlab = "control"
ylab <- "GATA6"
cutoff <- 2

all_data <- cbind(control_data, treatment_data)
## define the unchanged probes
No_change_num <- which((control_data-treatment_data) <cutoff |(control_data-treatment_data) > -cutoff)
No_change <- cbind(control_data[No_change_num], treatment_data[No_change_num])
nNo_change <- length(No_change_num) #umber ofunchanged probes

## define the methylated(up) probes
up_num <- which((control_data-treatment_data) <= -cutoff )
up <- cbind(control_data[up_num], treatment_data[up_num])
nup <- length(up_num) #umber of methylated probes

## define the demethylated(down) probes
down_num <- which((control_data-treatment_data) >= cutoff)
down <- cbind(control_data[down_num], treatment_data[down_num])
ndown <- length(down_num) #umber of demethylated probes

##setting of plot limit
data_range <- c(min(c(control_data, treatment_data)), max((c(control_data, treatment_data))))
plot_limit <- c(
    mean(data_range)-((abs(data_range[2]-data_range[1])/2)*1.1),
                mean(data_range) + ((abs(data_range[2] - data_range[1]) / 2) * 1.1)
)
##plotting
##plot all data as smooth scatter
smoothScatter(
    all_data[, 1],
    all_data[, 2],
    pch="",
    xlim = plot_limit,
    ylim = plot_limit,
    xlab = xlab,
    ylab = ylab,
    main=main,
    colramp = colorRampPalette(c("white", "gray18"))
)

par(new=T)
#methylated and demethylated probes
plot(0, 0, type = "n", xlim=plot_limit, ylim=plot_limit ,xlab = "", ylab = "")
points(up[,1], up[,2], pch=".", cex=4, col=rgb(0, 0.7, 0, alpha=0.2)) #methylated
points(down[, 1], down[, 2], pch = ".", cex = 4, col = rgb(0.8, 0, 0.2, alpha = 0.2)) #demethylated

#threshold lines
abline(-2,1,col="black",lty=2)
abline(2, 1, col = "black", lty = 2)

#numbers ofDMP
text(data_range[1], data_range[2], adj=0.3, paste("Methylated = ", nup, sep="") , col=rgb(0, 0.7, 0, alpha=0.9))
text(data_range[2], data_range[1], adj=0.8, paste("Demethylated = ",ndown,sep=""), col=rgb(0.8, 0, 0.2, alpha=0.9))

#===================================================================================================#

#Motf overrepresenation analysis (Fig. 3C)

#===================================================================================================#
g1 <- GATA6_motif_analysis(selDataMatrix, prefix = "GATA6_demethyl", MD = "Demethyl")
g2 <- GATA6_motif_analysis(selDataMatrix, prefix = "GATA6_methyl", MD = "Methyl")

g <- plot_grid(g2, g1, align =  "hv", axis = "tblr", nrow = 1)
print(g)