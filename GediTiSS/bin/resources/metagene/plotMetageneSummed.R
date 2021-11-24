library(ggplot2)

plotFile <- function(fileIn, fileOut) {
    df = read.delim(fileIn, header=TRUE)
    colorDf = data.frame(fillType=c(), xMin=c(), xMax=c())
    types = unique(df$Type)

    for (t in types) {
        indices=which(df$Type==t)
        colorDf = rbind(colorDf, data.frame(fillType=c(t), xMin=c(indices[1]), xMax=c(indices[length(indices)])))
    }

    gp = ggplot(data=df)
    gp = gp + geom_rect(data=colorDf, aes(ymin=-Inf, ymax=Inf, xmin=xMin, xmax=xMax, fill=fillType), alpha=0.4)
    gp = gp + geom_line(aes(x=Position, y=Value)) + ylim(0,yMaxValue)
    pdf(paste(fileOut, ".pdf", sep=""))
    print(gp)
    dev.off()
}

plotMerged <- function(df, colorDf) {
    gp = ggplot(data=df)
    gp = gp + geom_rect(data=colorDf, aes(ymin=-Inf, ymax=Inf, xmin=xMin, xmax=xMax, fill=fillType), alpha=0.4)
    gp = gp + geom_line(aes(x=Position, y=Value, color=groupi)) + ylim(0,yMaxValue)
    pdf(paste(prefixOut, "/combined.pdf", sep=""), width=10, height=6)
    print(gp)
    dev.off()
}

getColorDf = function(file) {
    df = read.delim(file, header=TRUE)
    colorDf = data.frame(fillType=c(), xMin=c(), xMax=c())
    types = unique(df$Type)

    for (t in types) {
        indices=which(df$Type==t)
        colorDf = rbind(colorDf, data.frame(fillType=c(t), xMin=c(indices[1]), xMax=c(indices[length(indices)]+1)))
    }

    return(colorDf)
}

colorDf = getColorDf(list.files(path=prefixIn, pattern="*.tsv$", full.names=TRUE)[1])
mergedDf = data.frame(Position=c(), Type=c(), Value=c(), groupi=c())

for (f in list.files(path=prefixIn, pattern="*.tsv$", full.names=TRUE)) {
    df = read.delim(f, header=TRUE)
    splitF = which(strsplit(f, "")[[1]]=="/")
    df$groupi = substr(f,splitF[length(splitF)]+1,nchar(f))
    mergedDf = rbind(mergedDf, df)
}
if (yMaxValue < 0) {
    yMaxValue = max(mergedDf$Value) * 1.2
}
plotMerged(mergedDf, colorDf)

for (f in list.files(path=prefixIn, pattern="*.tsv$", full.names=FALSE)) {
    plotFile(paste(prefixIn, "/", f, sep=""), paste(prefixOut, "/", f, sep=""))
}