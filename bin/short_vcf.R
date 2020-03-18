#!/usr/bin/env Rscript

DP4_to_freq<-function(DP4){
  if (!is.na(DP4)) {
    DP4_num=as.numeric(DP4)
    freq=(DP4_num[3]+DP4_num[4])/sum(DP4_num)*100
  } else
  {
    freq=NA
  }
  return(freq)
}

get_seg<-function(name_vector){
  SEG_split = name_vector[1:length(name_vector)-1]
  return(paste(SEG_split,collapse="_"))
}

suppressMessages(library('vcfR'))
library('readr')
args = commandArgs(trailingOnly=TRUE)
vcffile = args[1]

vcfdata <- read.vcfR(vcffile, verbose = FALSE)
vcfDP4 <- extract.info(vcfdata, element = 'DP4')
vcfvarscanfreq <- extract.gt(vcfdata, element = 'FREQ')[,"Sample1"]
vcfvarscanfreq_num <- as.numeric(sub("%", "",vcfvarscanfreq))
vcfDP4_split <- strsplit(vcfDP4,split=",")
vcfDP4freq <- round(suppressWarnings(sapply(vcfDP4_split,DP4_to_freq)),2)
vcfDP4freq[which(is.na(vcfDP4freq))] = vcfvarscanfreq_num[which(is.na(vcfDP4freq))]

frequency_table <- data.frame(getCHROM(vcfdata),getPOS(vcfdata),getREF(vcfdata),getALT(vcfdata),vcfDP4freq)
names(frequency_table)=c("SEGMENT","POSITION","REF","ALT","FREQUENCY")
cat(format_tsv(frequency_table))
