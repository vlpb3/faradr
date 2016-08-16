# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title:Build a function to read all the Index Statistics (samtools output) and combine them in one table for all samples


CombineCountTables <- function(path, output, fileSuffix) {
  listOfFiles = grep(fileSuffix, dir(path), value=TRUE)
  
  sampleNames <- strsplit(listOfFiles, "_", fixed=T) 
  
  df <- read.delim(paste(path, listOfFiles[1], sep="/"))
  colnames(df) <- c("ID", "Len", "Counts", "Nth")
  
  compileExper <- df[,c(1,3)]
  colnames(compileExper) <- c("ID", sampleNames[[1]][1])
  
  
  for (i in 2:length(listOfFiles)){
    df <- read.delim(paste(path, listOfFiles[i], sep="/"))
    colnames(df) <- c("ID", "Len", "Counts", "Nth")
    df <-  df[,c(1,3)]
    colnames(df) <- c("ID", sampleNames[[i]][1])
    compileExper <- merge(compileExper, df, by="ID")
     
  }
  
  write.table(compileExper, file=output, sep="\t", row.names = F,quote=F)
}

# CombineCountTables(path, fileSuffix ) 
