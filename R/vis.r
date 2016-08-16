#' Plot aligned reads over annotated region.
#' @param aln read alignments in GAlignments objet (GenomicTRanges)
#' @param annots genome annotations in GRanges object
#' @param plotfile path for output plot
#' @importFrom ggbio autoplot tracks
#' @importFrom ggplot2 aes
#' @export
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
#' @param annots genome annotations in GRanges object
#' @param plotfile path for output plot
#' @importFrom ggbio autoplot tracks
#' @importFrom ggplot2 aes
#' @importFrom IRanges coverage
#' @export
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
#' @importFrom ShortRead qa
#' @importFrom ggplot2 ggplot geom_boxplot 
PlotPerCycleQuality <- function(data.dir, file.pattern) {
    input.type <- "fastq"
    sreadq.qa <- qa(data.dir, file.pattern, type=input.type)
    perCycle <- sreadq.qa[["perCycle"]]
    perCycle.q <- perCycle$quality
    pcq <- perCycle.q[rep(seq_len(nrow(perCycle.q)), perCycle.q$Count), ]
    p <- ggplot(pcq, aes(factor(Cycle), Score))
    p <- p + geom_boxplot()
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

#' Plot number of reads in each sample.
#' 
#' Bar plot representing number of reads in each sample.
#' Samples are expected to be separate fastq files. 
#' @param fqc FastQA from package ShortRead
#' @return plot object
#' @importFrom ggplot2 ggplot geom_bar theme labs
#' @export
Aplot <- function(fqc) {
  df <- data.frame(nReads=fqc[['readCounts']]$read,
                   sample=row.names(fqc[['readCounts']]))
  p <- ggplot(df, aes(x=sample, y=nReads, fill=sample))
  p <- p + geom_bar(alpha=0.6, stat='identity')
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Number of Reads")
  return(p)
}
  
#' Plot number of reads in each sample for each grouping factor.
#' 
#' Bar plots representing number of reads in each sample.
#' There is one bar plot per grouping factor.
#' Samples are grouped and clored by grouping factor.
#' Grouping factors are extracted from design table.
#' Samples are expected to be separate fastq files. 
#' @param fqc FastQA from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return list of plot objects
#' @importFrom ggplot2 ggplot geom_bar theme labs
#' @importFrom stringr str_replace
#' @export
A.design.plot <- function(fqc, design.table) {
    # get names of groups from desing file
    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    df <- data.frame(nReads=fqc[['readCounts']]$read,
                     sampleid=str_replace(row.names(fqc[['readCounts']]), "\\.f(ast)?q", ""))

    df <- merge(df, design.table)
    plots <- lapply(groups, function(g.factor) {
                    # order by grouping factor
                    df <- df[order(df[g.factor]), ]
                    .e <- environment()
                    p <- ggplot(df, aes(x=sampleid, y=nReads, fill=factor(df[,g.factor])),
                                environment=.e)
                    p <- p + geom_bar(alpha=0.6, stat='identity')
                    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
                    p <- p + guides(fill=guide_legend(title=g.factor))
                    p <- p + labs(x="Sample", y="Number of Reads")
                    return(p)
                 })
    names(plots) <- groups
    return(plots)
}

#' Plot read lenght distribution 
#' 
#' Density plot representing freaquencies or readlenths.
#' Samples are explected to be separate fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @importFrom ggplot2 ggplot geom_density labs
#' @export
B1plot <- function(samples) {
  rlens <- lapply(samples, width)
  df <- data.frame(sample=rep(names(samples), lapply(rlens, length)),
                   rlen=unlist(rlens), row.names=NULL)
  q95rlen <- quantile(df$rlen, 0.95)
  p <- ggplot(df, aes(rlen, group=sample, colour=sample))
  p <- p + geom_density(alpha=I(0.4), adjust=3) + xlim(0, q95rlen)
  p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sample")
}

#' Plot read lenght distribution 
#' 
#' Density plot representing freaquencies or readlenths.
#' There is one plot per grouping factor.
#' Samples clored by grouping factor.
#' Grouping factors are extracted from design.table.
#' Samples are explected to be separate fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return list of plot objects
#' @importFrom ggplot2 ggplot geom_density labs
#' @importFrom stringr str_replace
#' @export
B1.design.plot <- function(samples, design.table) {

    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    rlens <- lapply(samples, width)
    df <- data.frame(sampleid=rep(names(samples), lapply(rlens, length)),
                     rlen=unlist(rlens), row.names=NULL)
    df <- merge(df, design.table)
    q95rlen <- quantile(df$rlen, 0.95)

    plots <- lapply(groups, function(g.factor){
           .e <- environment()
           p <- ggplot(df, aes(rlen, group=sampleid, colour=factor(df[,g.factor])),
                       environment=.e)
           p <- p + geom_density(alpha=I(0.4), adjust=3) + xlim(0, q95rlen)
           p <- p + guides(colour=guide_legend(title=g.factor)) 
           p <- p + labs(x="Read Length", y="Fraction of Reads")
                     })
    names(plots) <- groups
    return(plots)
}

#' Plot fraction of reads with particular lengh or longer.
#' 
#' Line plot showing fraction of reads on y axis
#' and minimal read length on x axis.
#' One line per sample.
#' Samples are explected to be separage fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @importFrom ggplot2 ggplot geom_line xlim labs
#' @import dplyr
#' @export
B2plot <- function(samples) {
    rlens <- lapply(samples, width)
    df <- data.frame(sampleid=rep(names(rlens), lapply(rlens, length)),
                     rlen=unlist(rlens), row.names=NULL)
    q95rlen <- quantile(df$rlen, 0.95)

    df <- df %>% group_by(sampleid, rlen) %>% summarise(count=n())
    df <- df %>% group_by(sampleid) %>% mutate(cum.count=cumsum(count),
                                               frac.count=cum.count/sum(count))

    p <- ggplot(df, aes(group=sampleid, colour=sampleid ))
    p <- p + geom_line(aes(rlen, 1-frac.count), alpha=0.4)
    p <- p + xlim(min(df$rlen), q95rlen)
    p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sampleid")
}

#' Plot fraction of reads with particular lengh or longer.
#' 
#' Line plot showing fraction of reads on y axis
#' and minimal read length on x axis.
#' One line per sample.
#' There is one plot per grouping factor.
#' Samples clored by grouping factor.
#' Grouping factors are extracted from design.table.
#' Samples are explected to be separage fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return plots list of plot objects
#' @importFrom ggplot2 ggplot geom_line xlim labs
#' @import dplyr
#' @export
B2.design.plot <- function(samples, design.table) {
    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    rlens <- lapply(samples, width)
    df <- data.frame(sampleid=rep(names(rlens), lapply(rlens, length)),
                     rlen=unlist(rlens), row.names=NULL)
    q95rlen <- quantile(df$rlen, 0.95)
    df <- df %>% group_by(sampleid, rlen) %>% summarise(count=n())
    df <- df %>% group_by(sampleid) %>% mutate(cum.count=cumsum(count),
                                               frac.count=cum.count/sum(count))
    df <- merge(df,design.table)

    plots <- lapply(groups, function(g.factor){
        .e <- environment()
        p <- ggplot(df, aes(group=sampleid, colour=factor(df[ ,g.factor])),
                    environment=.e)
        p <- p + geom_line(aes(rlen, 1-frac.count), alpha=0.4)
        p <- p + xlim(min(df$rlen), q95rlen)
        p <- p + guides(colour=guide_legend(title=g.factor)) 
        p <- p + labs(x="Read Length", y="Fraction of Reads", colour="Sampleid")
    })
    names(plots) <- groups
    return(plots)
}

#' Plot mean read quality distribution per sample.
#' 
#' Boxplot showing the distribution of mean read quality.
#' One boxplot per sample.
#' Samples are explected to be separage fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @importFrom ggplot2 ggplot geom_boxplot theme labs
#' @export
C1plot <- function(samples) {
  calc.qmeans <- function(fq) {
    qm <- as(quality(fq), "matrix")
    row.means <- rowMeans(qm, na.rm=T)
    return(row.means)
  }
  qmeans <- lapply(samples, calc.qmeans)
  df <- data.frame(sample=rep(names(qmeans), lapply(qmeans, length)),
                   qmeans=unlist(qmeans), row.names=NULL)
  p <- ggplot(df, aes(factor(sample), qmeans, fill=sample))
  p <- p + geom_boxplot(alpha=0.6, outlier.size=0)
  p <- p + theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  p <- p + labs(x="Sample", y="Mean Read Quality")
}

#' Plot mean read quality distribution per sample.
#' 
#' Boxplot showing the distribution of mean read quality.
#' One boxplot per sample.
#' There is one plot per grouping factor.
#' Samples clored by grouping factor.
#' Grouping factors are extracted from design.table.
#' Samples are explected to be separage fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return list of plot objects
#' @importFrom ggplot2 ggplot geom_boxplot theme labs
#' @export
C1.design.plot <- function(samples, design.table) {

    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    calc.qmeans <- function(fq) {
        qm <- as(quality(fq), "matrix")
        row.means <- rowMeans(qm, na.rm=T)
        return(row.means)
    }
    qmeans <- lapply(samples, calc.qmeans)
    df <- data.frame(sampleid=rep(names(qmeans), lapply(qmeans, length)),
                     qmeans=unlist(qmeans), row.names=NULL)
    df <- merge(df, design.table)
    plots <- lapply(groups, function(g.factor) {
                    df <- df[order(df[g.factor]), ]
                    .e <- environment()
                    p <- ggplot(df, aes(factor(sampleid), qmeans, fill=factor(df[ ,g.factor])),
                                environment=.e)
                    p <- p + geom_boxplot(alpha=0.6, outlier.size=0)
                    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
                    p <- p + guides(fill=guide_legend(title=g.factor))
                    p <- p + labs(x="Sampleid", y="Mean Read Quality")
                     })
    names(plots) <- groups
    return(plots)
}

#' Plot mean base quality at particualar position in the read.
#' 
#' Line plot showing mean quality of the bases per position in the read.
#' One line per sample.
#' Samples are explected to be separage fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @import dplyr
#' @importFrom ggplot2 ggplot geom_line xlim labs
#' @export
C2plot <- function(samples) {
    calc.qmeans <- function(fq) {
        qm <- as(quality(fq), "matrix")
        col.means <- colMeans(qm, na.rm=T)
    }
    q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
    qmeans <- lapply(samples, calc.qmeans)
    df <- data.frame(sampleid=rep(names(qmeans), lapply(qmeans, length)),
                     qmeans=unlist(qmeans), row.names=NULL)

    df <- df %>% group_by(sampleid) %>% mutate(ones=1, pos=cumsum(ones))

    p <- ggplot(df, aes(pos, qmeans, group=sampleid, colour=sampleid))
    p <- p + geom_line(alpha=0.4) + xlim(0, q95rlen)
    p <- p + labs(x="Position in the Read", y="Mean Base Quality", colour="Sampleid")
}

#' Plot mean base quality at particualar position in the read.
#' 
#' Line plot showing mean quality of the bases per position in the read.
#' One line per sample.
#' Samples are explected to be separage fastq files. 
#' There is one plot per grouping factor.
#' Samples colored by grouping factor.
#' Grouping factors are extracted from design.table.
#' @param samples ShortReadQ object from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return list of plot objects
#' @import dplyr
#' @importFrom ggplot2 ggplot geom_line xlim labs
#' @export
C2.design.plot <- function(samples, design.table) {
    
    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    calc.qmeans <- function(fq) {
        qm <- as(quality(fq), "matrix")
        col.means <- colMeans(qm, na.rm=T)
    }
    q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
    qmeans <- lapply(samples, calc.qmeans)
    df <- data.frame(sampleid=rep(names(qmeans), lapply(qmeans, length)),
                     qmeans=unlist(qmeans), row.names=NULL)
    df <- df %>% group_by(sampleid) %>% mutate(ones=1, pos=cumsum(ones))
    df <- merge(df, design.table)
    plots <- lapply(groups, function(g.factor){
        .e <- environment()
        p <- ggplot(df, aes(pos, qmeans, group=sampleid, colour=factor(df[ ,g.factor])),
                    environment=.e)
        p <- p + geom_line(alpha=0.4) + xlim(0, q95rlen)
        p <- p + guides(colour=guide_legend(title=g.factor)) 
        p <- p + labs(x="Position in the Read", y="Mean Base Quality")
                     }) 
    names(plots) <- groups
    return(plots)
}

#' Plot fraction of reads with partucular quality or higher per position.
#'
#' 4 line plots, for qualities 18, 20, 24 and 28, showing fraction 
#' of reads with one of those qualities per position in the read.
#' One line per sample.
#' Samples are explected to be separage fastq files. 
#' There is one plot per grouping factor.
#' Samples colored by grouping factor.
#' Grouping factors are extracted from design.table.
#' @param samples ShortReadQ object from package ShortRead
#' @param fqc FastQA from package ShortRead
#' @param design.table data.frame holds information about experimantal design 
#' @return list of plot objects
#' @import dplyr
#' @importFrom ggplot2 ggplot geom_line facet_wrap ylim labs
#' @export
C3.design.plot <- function(samples, fqc, design.table) {
    
    groups <- names(design.table)
    groups <- groups[groups != "sampleid"]

    border.quals <- c(18, 20, 24, 28)
    q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
    pcq <- fqc[["perCycle"]]$quality
    pcq$sampleid <- str_replace(pcq$lane, "\\.f(ast)?q", "")
    pcq <-  merge(pcq, design.table, by='sampleid')
    pcq <- pcq %>% group_by(sampleid, Cycle) %>% mutate(CycleCounts=sum(Count),
                                                        CountFrac=Count/CycleCounts,
                                                        ScoreCumSum=cumsum(CountFrac)) %>% ungroup()
    pcq <- pcq[pcq$Cycle <= q95rlen, ]
    subpcq <- pcq[pcq$Score %in% border.quals, ]
    plots <- lapply(c(groups), function(g.factor){
        .e <- environment()
        p <- ggplot(subpcq, environment=.e)
        p <- p + geom_line(aes_string(x='Cycle', y='1 - ScoreCumSum',
                                      group='sampleid', colour=g.factor),
                           alpha=I(0.4))
        p <- p + facet_wrap(~Score) + ylim(0,1)
        p <- p + guides(colour=guide_legend(title=g.factor)) 
        p <- p + labs(x="Position in the Read", y="Fraction of Reads")
                 })
    names(plots) <- groups
    return(plots)
}

#' Plot fraction of reads with partucular quality or higher per position.
#'
#' 4 line plots, for qualities 18, 20, 24 and 28, showing fraction 
#' of reads with one of those qualities per position in the read.
#' One line per sample.
#' Samples are explected to be separage fastq files. 
#' @param fqc FastQA from package ShortRead
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @import dplyr
#' @importFrom ggplot2 ggplot geom_line facet_wrap ylim labs
#' @export
C3plot <- function(fqc, samples) {
    border.quals <- c(18, 20, 24, 28)
    q95rlen <- quantile(unlist(sapply(samples, width)), 0.95)
    pcq <- fqc[["perCycle"]]$quality
    pcq$sampleid <- str_replace(pcq$lane, "\\.f(ast)?q", "")

    pcq <- pcq %>% group_by(sampleid, Cycle) %>% mutate(CycleCounts=sum(Count),
                                                        CountFrac=Count/CycleCounts,
                                                        ScoreCumSum=cumsum(CountFrac)) %>% ungroup()
    pcq <- pcq[pcq$Cycle <= q95rlen, ]
    subpcq <- pcq[pcq$Score %in% border.quals, ]
    p <- ggplot(subpcq)
    p <- p + geom_line(aes(Cycle, 1-ScoreCumSum, group=sampleid, colour=sampleid), alpha=I(0.4))
    p <- p + facet_wrap(~Score) + ylim(0,1)
    p <- p + labs(x="Position in the Read", y="Fraction of Reads", colour="Sampleid")
}

#' Plot base frequnecy for first 30 nt.
#'
#' Linepolots representing frequencies of bases, one line per base,
#' for first 30 nt in the read (from 5')
#' One plot per sample.
#' Each sample is a fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line scale_x_discrete facet_wrap labs
#' @importFrom ggplot2 element_rect
#' @import dplyr
#' @export
D1plot <- function(samples) {
    countLetterFreq <- function(fq) {
        bases <- c("A", "C", "G", "T")
        sr <- sread(fq)
        alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
        alpha.freq <- data.frame(alpha.freq)
        colnames(alpha.freq) <- 1:30
        alpha.freq$base <- rownames(alpha.freq)
        df <- melt(alpha.freq, id=c("base"))
        colnames(df) <- c("base", "pos", "count")
        df <- df %>% group_by(pos) %>% mutate(total.count=sum(count),
                                              frac.count=count/total.count)
        # df <- ddply(df, .(pos), mutate,
        #             total.count=sum(count),
        #             frac.count=count/total.count)
        return(df)
    }
    df.list <- lapply(samples, countLetterFreq)
    df <- do.call("rbind", df.list)
    df$sampleid <- rep(names(samples),each=4*30)
#    df$sampleid <-sapply(rownames(df),
#                         function(rn) {str_split(rn, "\\.")[[1]][1]},
#                         simplify=TRUE, USE.NAMES=FALSE) 
     
    p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
    p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(0, 30, 5))
    p <- p + ylim(0,1)
    p <- p + facet_wrap(~sampleid)
    p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
    p <- p + theme(panel.background=element_rect(fill="white", colour="grey"))
}

#' Plot base frequnecy for last 30 nt.
#'
#' Lineplots representing frequencies of bases, one line per base,
#' for last 30 nt in the read (from 3')
#' One plot per sample.
#' Each sample is a fastq files. 
#' @param samples ShortReadQ object from package ShortRead
#' @return plot object
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line scale_x_discrete facet_wrap labs
#' @importFrom ggplot2 element_rect
#' @import dplyr
#' @export
D2plot <- function(samples) {
    countLetterFreq <- function(fq) {
        sr <- reverse(sread(fq))
        bases <- c("A", "C", "G", "T")
        alpha.freq <- alphabetByCycle(sr, bases)[, 1:30]
        alpha.freq <- data.frame(alpha.freq)
        colnames(alpha.freq) <- -1:-30
        alpha.freq$base <- rownames(alpha.freq)
        df <- melt(alpha.freq, id=c("base"))
        colnames(df) <- c("base", "pos", "count")
        df <- df %>% group_by(pos) %>% mutate(total.count=sum(count),
                                              frac.count=count/total.count)
        return(df)
    }
    
    df.list <- lapply(samples, countLetterFreq)
    df <- do.call("rbind", df.list)
    df$sampleid <- rep(names(samples),each=4*30)
#    df$sampleid <-sapply(rownames(df),
#                         function(rn) {str_split(rn, "\\.")[[1]][1]},
#                         simplify=TRUE, USE.NAMES=FALSE) 

  p <- ggplot(df, aes(pos, frac.count, group=base, colour=base))
  p <- p + geom_line(alpha=0.4) + scale_x_discrete(breaks=seq(-30, 0, 5))
  p <- p + ylim(0,1)
  p <- p + facet_wrap(~sampleid)
  p <- p + labs(x="Position in the Read", y="Base Frequency", colour="Base")
  p <- p + theme(panel.background=element_rect(fill="white", colour="white"))
}
