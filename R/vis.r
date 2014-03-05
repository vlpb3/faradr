library(ShortRead)
library(ggbio)
library(ggplot2)
library(plyr)

#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annot genome annotations in GRanges object
#' @param plotfile path for output plot
PlotAlignAnnot <- function(aln, annots, plotfile) {
    reads.track <- autoplot(aln)
    annots.track <- autoplot(annots,
                             aes(color=group, fill=group))
    png(filename=plotfile,
        width=1024, height=512)  
    print(tracks(alignments=reads.track,
                 annotation=annots.track,
                 heights=c(5,1)))  
    dev.off()
}

#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annot genome annotations in GRanges object
#' @param plotfile path for output plot
PlotCoverageAnnot <- function(aln, annots, plotfile) {
    cov <- coverage(aln)
    coverage.track <- autoplot(cov, binwidth=10)
    
    annots.track <- autoplot(annots,
                             aes(color=group, fill=group))
    png(filename=plotfile,
        width=1024, height=512)  
    print(tracks(coverage=coverage.track,
                 annotation=annots.track,
                 heights=c(5,1)))  
    dev.off()
}

#' Plot quality scores per cycle
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @return plot object
PlotPerCycleQuality <- function(data.dir, file.pattern) {
    require(ggplot2)
    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    perCycle.q <- perCycle$quality
    pcq <- perCycle.q[rep(seq_len(nrow(perCycle.q)), perCycle.q$Count), ]
    p <- ggplot(pcq, aes(factor(Cycle), Score))
    p <- p + geom_boxplot()
    return(p)
}

#' Plot cumulative qualities 
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @return plot 
PlotPerCycleCumQuality <- function(data.dir, file.pattern) {
    require(ggplot2)
    require(plyr)
    fastq.files <- list.files(data.dir, file.pattern)
    q95rlen.df <- adply(fastq.files, 1, function(f, data.dir) {
      fq <- readFastq(data.dir, f)
      q95rlen <- quantile(width(fq), .95)
      data.frame(lane=c(f), q95rlen=c(q95rlen))
    }, data.dir=data.dir)
    q95rlen.df <- subset(q95rlen.df, select = c(lane, q95rlen))

    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    pcq <- perCycle$quality
    pcq <- ddply(pcq, .(lane, Cycle), mutate, CycleCounts=sum(Count))    
    pcq$CountFrac <- pcq$Count / pcq$CycleCounts
    pcq <- ddply(pcq, .(lane, Cycle), mutate, ScoreCumSum=cumsum(CountFrac))
    pcq <- ddply(pcq, .(lane), mutate, LaneCounts=CycleCounts[1])
    pcq <- merge(x=pcq, y=q95rlen.df, by="lane", all=TRUE)
    pcq <- pcq[pcq$Cycle <= pcq$q95rlen, ]
    subpcq <- pcq[pcq$Score %in% c(15, 20, 25, 30), ]
    p <- ggplot(subpcq)
    p <- p + geom_line(aes(Cycle, 1-ScoreCumSum, group=Score, colour=factor(Score)), alpha=I(0.8))
    p <- p + geom_area(aes(Cycle, CycleCounts / LaneCounts), position="identity", alpha=I(0.1))
    p <- p + facet_wrap(~lane, scales="free_x")
    return(p)
}

#' Plot mean quality per read
#' @param data.dir path to data dir
#' @param file.pattern pattern of input fastq
#' @return plot object
PlotMeanReadQuality <- function(data.dir, file.pattern) {
  count.qmeans <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    qm <- as(quality(fq), "matrix")
    row.means <- rowMeans(qm, na.rm=T)
    qmeans <- data.frame(mean.qual=row.means, lane=rep(s, length(row.means)))
    return(qmeans)
  }
  samples <- list.files(data.dir, file.pattern)
  qmeans <- adply(samples, 1, count.qmeans, data.dir=data.dir)
  p <- ggplot(qmeans, aes(x=mean.qual))
  p <- p + geom_histogram(aes(y=..density..),
                          alpha=.1,
                          fill="green",
                          colour="darkgreen",
                          binwidth=.5)
  p <- p + geom_density() + facet_wrap(~lane)
  return(p)
}

#' Plot read length distribution
#' @param data.dir path to data dir
#' @param file.pattern pattern of input fastq
#' @return plot object 
PlotReadLengthDistribution <- function(data.dir, file.pattern) {
  samples <- list.files(data.dir, file.pattern)
  count.rlens <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    rlens <- width(fq)
    rlen95 <- quantile(width(fq), .95)
    rlens <- rlens[rlens <= rlen95]
    rlens.df <- data.frame(read.len=rlens, lane=rep(s, length(rlens)))
    return(rlens.df)
  }
  
  rlens <- adply(samples, 1, count.rlens, data.dir=data.dir)
  p <- ggplot(rlens, aes(x=read.len))
  p <- p + geom_histogram(aes(y=..density..),
                          alpha=0.02,
                          fill="green",
                          colour="darkgreen",
                          binwidth=1)
  p <- p + geom_density(adjust=3) + facet_wrap(~lane)
  return(p)
}

#' Plot base call per cycle
#' @param data.dir string path to data
#' @param file.pattern string pattern for input files
#' @param type stirng input type, 'fastq' is default value
#' @return plot 
PlotPerCycleBaseCalls <- function(data.dir, file.pattern, type="fastq") {
    require(ggplot2)  
    require(plyr)
    
    sreadq.qa <- qa(data.dir, file.pattern, type=type)
    perCycleBaseCall <- sreadq.qa[["perCycle"]]$baseCall
    perCycleCounts  <- ddply(perCycleBaseCall, c("Cycle"), summarise,
                           countsInCycle = sum(Count))
    perCycleBaseCall$callsFraction <- perCycleBaseCall$Count / perCycleCounts$countsInCycle[perCycleBaseCall$Cycle] 
    p <- ggplot(perCycleBaseCall, aes(Cycle, callsFraction))
    p <- p + geom_point(aes(colour=factor(Base)), alpha=0.7)
    p <- p + geom_line(aes(colour=factor(Base)), alpha=0.7) 
    return(p)
}

Aplot <- function(data.dir, file.pattern) {
  samples <- list.files(data.dir, file.pattern)    
  countReads <- function(s, data.dir) {
    fq <- readFastq(data.dir, s)  
    return(length(fq))
  }
  counts <- sapply(samples, countReads, data.dir=data.dir)
  df <- data.frame(nReads=counts, sample=names(counts))
  p <- ggplot(df, aes(x=sample, y=nReads, fill=sample))
  p <- p + geom_bar(alpha=0.6, stat='identity')
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Number of Reads")
  return(p)
}

Aplot2 <- function(fqc) {
  df <- data.frame(nReads=fqc[['readCounts']]$read,
                   sample=row.names(fqc[['readCounts']]))
  p <- ggplot(df, aes(x=sample, y=nReads, fill=sample))
  p <- p + geom_bar(alpha=0.6, stat='identity')
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Number of Reads")
  return(p)
}
  
B1plot <- function(data.dir, file.pattern){
  samples <- list.files(data.dir, file.pattern)
  get.rlens <- function(s, data.dir) {
    sample.size = 1000000
    fq <- yield(FastqSampler(file.path(data.dir, s), n=sample.size))  
    rlens <- width(fq)
    rlens.df <- data.frame(read.length=rlens, sample=rep(s, length(rlens)))
    return(rlens.df)
  }
  df <- adply(samples, 1, get.rlens, data.dir=data.dir)
  rlen.limit <- quantile(df$read.length, 0.95)
  p <- ggplot(df, aes(read.length, group=sample, colour=sample))
  p <- p + geom_density(alpha=I(0.4), adjust=3) + xlim(0, rlen.limit)
  p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sample")
}

B1plot2 <- function(samples) {
  rlens <- lapply(samples, width)
  df <- data.frame(sample=rep(names(samples), lapply(rlens, length)),
                   rlen=unlist(rlens), row.names=NULL)
  q95rlen <- quantile(df$rlen, 0.95)
  p <- ggplot(df, aes(rlen, group=sample, colour=sample))
  p <- p + geom_density(alpha=I(0.4), adjust=3) + xlim(0, q95rlen)
  p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sample")
}

B2plot <- function(data.dir, file.pattern) {
  samples <- list.files(data.dir, file.pattern)
  get.rlens <- function(s, data.dir) {
    sample.size = 1000000
    fq <- yield(FastqSampler(file.path(data.dir, s), n=sample.size))  
    rlens <- width(fq)
    rlens.df <- data.frame(read.length=rlens, sample=rep(s, length(rlens)))
    return(rlens.df)
  }
  df <- adply(samples, 1, get.rlens, data.dir=data.dir)
  # count read length accurences per sample
  df.sum <- ddply(df, .(sample),
                  summarise,
                  rlen=as.integer(names(table(read.length))),
                  count=table(read.length))
  # cunt counts in sample (some samples have less then 1M reads)
  df.sum <- ddply(df.sum, .(sample), mutate,
                  totall=sum(count),
                  frac=count/totall,
                  cums=cumsum(frac))  
  
  rlen.limit <- quantile(df$read.length, 0.95)
  df.sum <- df.sum[df.sum$rlen <= rlen.limit, ]
  p <- ggplot(df.sum, aes(group=sample, colour=sample ))
  p <- p + geom_line(aes(rlen, 1-cums), alpha=0.4)
  p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sample")
  return(p)
}

B2plot2 <- function(samples) {
  rlens <- lapply(samples, width)
  df <- data.frame(sample=rep(names(rlens), lapply(rlens, length)),
                   rlen=unlist(rlens), row.names=NULL)
  q95rlen <- quantile(df$rlen, 0.95)
  df <- ddply(df, .(sample, rlen), transform, counts=length(rlen))
  df <- df[!duplicated(df), ]
  df <- ddply(df, .(sample), mutate,
              cum.counts=cumsum(counts),
              frac.counts=cum.counts/sum(counts))
  p <- ggplot(df, aes(group=sample, colour=sample ))
  p <- p + geom_line(aes(rlen, 1-frac.counts), alpha=0.4)
  p <- p + xlim(min(df$rlen), q95rlen)
  p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sample")
}

C1plot <- function(data.dir, file.pattern) {
  count.qmeans <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    qm <- as(quality(fq), "matrix")
    row.means <- rowMeans(qm, na.rm=T)
    qmeans <- data.frame(mean.qual=row.means, sample=rep(s, length(row.means)))
    return(qmeans)
  }
  samples <- list.files(data.dir, file.pattern)
  qmeans <- adply(samples, 1, count.qmeans, data.dir=data.dir)
  p <- ggplot(qmeans, aes(factor(sample), mean.qual, fill=sample))
  p <- p + geom_boxplot(alpha=0.6)
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Mean Read Quality")
  return(p)
}

C1plot2 <- function(samples) {
  calc.qmeans <- function(fq) {
    qm <- as(quality(fq), "matrix")
    row.means <- rowMeans(qm, na.rm=T)
    return(row.means)
  }
  qmeans <- lapply(samples, calc.qmeans)
  df <- data.frame(sample=rep(names(qmeans), lapply(qmeans, length)),
                   qmeans=unlist(qmeans), row.names=NULL)
  p <- ggplot(df, aes(factor(sample), qmeans, fill=sample))
  p <- p + geom_boxplot(alpha=0.6)
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Mean Read Quality")
}

C2plot <- function(data.dir, file.pattern) {
  count.qmeans <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    q95rlen <- quantile(width(fq), 0.95)
    qm <- as(quality(fq), "matrix")
    col.means <- colMeans(qm, na.rm=T)
    qmeans <- data.frame(pos=1:length(col.means),
                         mean.qual=col.means,
                         sample=rep(s, length(col.means)),
                         q95rlen=rep(q95rlen, length(col.means)))
    return(qmeans)
  }
  samples <- list.files(data.dir, file.pattern)
  qmeans <- adply(samples, 1, count.qmeans, data.dir=data.dir)
  rlen.limit <- max(qmeans$q95rlen)  
  p <- ggplot(qmeans, aes(pos, mean.qual, group=sample, colour=sample))
  p <- p + geom_line(alpha=0.4) + xlim(0,rlen.limit)
  p <- p + labs(x="Position in the Read", y="Mean Base Quality", colour="Sample")
  return(p)
}

C2plot2 <- function(samples) {
  calc.qmeans <- function(fq) {
    qm <- as(quality(fq), "matrix")
    col.means <- colMeans(qm, na.rm=T)
  }
  q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
  qmeans <- lapply(samples, calc.qmeans)
  df <- data.frame(sample=rep(names(qmeans), lapply(qmeans, length)),
                   qmeans=unlist(qmeans), row.names=NULL)
  df <- ddply(df, .(sample), transform, pos=1:length(qmeans))  
  p <- ggplot(df, aes(pos, qmeans, group=sample, colour=sample))
  p <- p + geom_line(alpha=0.4) + xlim(0, q95rlen)
  p <- p + labs(x="Position in the Read", y="Mean Base Quality", colour="Sample")
}

C3plot <- function(data.dir, file.pattern) {
    fastq.files <- list.files(data.dir, file.pattern)
    q95rlen.df <- adply(fastq.files, 1, function(f, data.dir) {
      fq <- readFastq(data.dir, f)
      q95rlen <- quantile(width(fq), .95)
      data.frame(lane=c(f), q95rlen=c(q95rlen))
    }, data.dir=data.dir)
    q95rlen.df <- subset(q95rlen.df, select = c(lane, q95rlen))

    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    pcq <- perCycle$quality
    pcq <- ddply(pcq, .(lane, Cycle), mutate, CycleCounts=sum(Count))    
    pcq$CountFrac <- pcq$Count / pcq$CycleCounts
    pcq <- ddply(pcq, .(lane, Cycle), mutate, ScoreCumSum=cumsum(CountFrac))
    pcq <- ddply(pcq, .(lane), mutate, LaneCounts=CycleCounts[1])
    pcq <- merge(x=pcq, y=q95rlen.df, by="lane", all=TRUE)
    pcq <- pcq[pcq$Cycle <= pcq$q95rlen, ]
    subpcq <- pcq[pcq$Score %in% c(18, 20, 24, 28), ]
    p <- ggplot(subpcq)
    p <- p + geom_line(aes(Cycle, 1-ScoreCumSum, group=lane, colour=lane), alpha=I(0.4))
    p <- p + facet_wrap(~Score) + ylim(0,1)
    p <- p + labs(x="Position in the Read", y="Fraction of Reads", colour="Sample")
}

C3plot2 <- function(fqc, samples) {
  q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
  pcq <- fqc[["perCycle"]]$quality
  pcq <- ddply(pcq, .(lane, Cycle), mutate,
               CycleCounts=sum(Count),
               CountFrac=Count/CycleCounts,
               ScoreCumSum=cumsum(CountFrac))    
  pcq <- ddply(pcq, .(lane), transform, LaneCounts=CycleCounts[1])
  pcq <- pcq[pcq$Cycle <= q95rlen, ]
  subpcq <- pcq[pcq$Score %in% c(18, 20, 24, 28), ]
  p <- ggplot(subpcq)
  p <- p + geom_line(aes(Cycle, 1-ScoreCumSum, group=lane, colour=lane), alpha=I(0.4))
  p <- p + facet_wrap(~Score) + ylim(0,1)
  p <- p + labs(x="Position in the Read", y="Fraction of Reads", colour="Sample")
}

D1plot <- function(data.dir, file.pattern) {
  require(reshape)
  countLetterFreq <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    sr <- sread(fq)
    bases <- c("A", "C", "G", "T")
    alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
    alpha.freq <- data.frame(alpha.freq)
    colnames(alpha.freq) <- 1:30
    alpha.freq$base <- rownames(alpha.freq)
    alpha.freq$lane <- rep(s, 4)
    df <- melt(alpha.freq, id=c("base", "lane"))
    colnames(df) <- c("base", "lane", "pos", "count")
    df <- ddply(df, .(lane, pos), mutate,
                               total.count=sum(count),
                               frac.count=count/total.count)
    return(df)
  }
  samples <- list.files(data.dir, file.pattern)
  df <- adply(samples, 1, countLetterFreq, data.dir=data.dir)
  p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
  p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(0, 30, 5))
  p <- p + facet_wrap(~lane)
  p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
}

D1plot2 <- function(samples) {
  require(reshape)
  countLetterFreq <- function(fq) {
    sr <- sread(fq)
    bases <- c("A", "C", "G", "T")
    alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
    alpha.freq <- data.frame(alpha.freq)
    colnames(alpha.freq) <- 1:30
    alpha.freq$base <- rownames(alpha.freq)
    df <- melt(alpha.freq, id=c("base"))
    colnames(df) <- c("base", "pos", "count")
    df <- ddply(df, .(pos), mutate,
                               total.count=sum(count),
                               frac.count=count/total.count)
    return(df)
  }
  df <- adply(samples, 1, countLetterFreq)
  sample.labels <- as.vector(sapply(names(samples), function(l) rep(l, 3*40)))
  df$sample <- sample.labels
  p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
  p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(0, 30, 5))
  p <- p + facet_wrap(~sample)
  p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
}

D2plot <- function(data.dir, file.pattern) {
  require(reshape)
  countLetterFreq <- function(s, data.dir) {
    fq <- yield(FastqSampler(file.path(data.dir, s), n=1000000))  
    sr <- reverse(sread(fq))
    bases <- c("A", "C", "G", "T")
    alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
    alpha.freq <- data.frame(alpha.freq)
    colnames(alpha.freq) <- -1:-30
    alpha.freq$base <- rownames(alpha.freq)
    alpha.freq$lane <- rep(s, 4)
    df <- melt(alpha.freq, id=c("base", "lane"))
    colnames(df) <- c("base", "lane", "pos", "count")
    df <- ddply(df, .(lane, pos), mutate,
                               total.count=sum(count),
                               frac.count=count/total.count)
    return(df)
  }
  samples <- list.files(data.dir, file.pattern)
  df <- adply(samples, 1, countLetterFreq, data.dir=data.dir)
  p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
  p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(-30, 0, 5))
  p <- p + facet_wrap(~lane)
  p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
}

D2plot2 <- function(samples) {
  require(reshape)
  countLetterFreq <- function(fq) {
    sr <- reverse(sread(fq))
    bases <- c("A", "C", "G", "T")
    alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
    alpha.freq <- data.frame(alpha.freq)
    colnames(alpha.freq) <- -1:-30
    alpha.freq$base <- rownames(alpha.freq)
    df <- melt(alpha.freq, id=c("base"))
    colnames(df) <- c("base", "pos", "count")
    df <- ddply(df, .(pos), mutate,
                               total.count=sum(count),
                               frac.count=count/total.count)
    return(df)
  }
  df <- adply(samples, 1, countLetterFreq)
  sample.labels <- as.vector(sapply(names(samples), function(l) rep(l, 3*40)))
  df$sample <- sample.labels
  p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
  p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(-30, 0, 5))
  p <- p + facet_wrap(~sample)
  p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
}