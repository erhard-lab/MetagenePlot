library(ComplexHeatmap)
library(circlize)
library(viridis)

plotFile <- function(fileIn, fileOut, binNames, binSizes, names, noHmNorm) {
    df = read.delim(fileIn, header=FALSE)

    colSplitTmp = c()
    for (i in c(1:length(binNames))) {
        colSplitTmp = c(colSplitTmp, rep(binNames[i], binSizes[i]))
    }
    colSplit = factor(colSplitTmp, levels=binNames)

    mat = as.matrix(df)
#     mat = t(apply(mat, 1, function(x) x/max(x)))
    col = colorRamp2(breaks=c(0, seq(1E-8, 1, length.out=9)), colors=c("lightgrey", rev(magma(9))))

    colnames(mat) = rep("", dim(mat)[2])

    hm = Heatmap(name="Coverage",mat,cluster_rows=F,cluster_columns=F,col=col,column_split=colSplit,column_gap=unit(c(1),"mm"))
    pdf(paste(fileOut, ".pdf", sep=""))
    if (length(names) != 0) {
        draw(hm + rowAnnotation(rn = anno_text(names, gp = gpar(fontsize = 5))))
    } else {
        draw(hm)
    }
    dev.off()
}

for (f in list.files(path=prefixIn, pattern="*.tsv$", full.names=FALSE)) {
    plotFile(paste(prefixIn, "/", f, sep=""), paste(prefixOut, "/", f, sep=""), binNames, binSizes, names, noHmNorm)
}
