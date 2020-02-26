#!/usr/bin/env Rscript

get_ref<-function(csvfile){
  return (csvfile["Segment",1])
}

extract_reference<-function(argref,argcsvfile,argcsvrnames){
  noname_frame = bind_cols(argcsvfile[which(sapply(argcsvfile,get_ref)==argref)])
  namedframe = cbind(argcsvrnames,noname_frame)
  return(namedframe)
}

add_reference_segment<-function(covfile,filename){
  splitfilename=unlist(strsplit(unlist(strsplit(filename,'[.]'))[2],"_"))
  sampleref=splitfilename[1]
  sampleseg=splitfilename[2]
  tempframe = data.frame(as.character(sampleseg),stringsAsFactors = FALSE)
  colnames(tempframe) = as.character(sampleref)
  rownames(tempframe) = "Segment"
  remframe = covfile[(row.names(covfile)!="Depth Coverage Level")&(row.names(covfile)!="Zero Coverage"),,drop=FALSE]
  colnames(remframe) = sampleref
  return(rbind.data.frame(tempframe,remframe,stringsAsFactors = FALSE))
}

library(dplyr)

args = commandArgs(trailingOnly=TRUE)
sampleid=args[1]
covsumext=args[2]

file_list<-list.files(pattern=paste(sampleid,"*",covsumext,sep="."))
csvlist<- lapply(file_list,read.csv,sep="\t", row.names = 1, stringsAsFactors = FALSE)
csvproc<-mapply(FUN = add_reference_segment, csvlist, file_list, SIMPLIFY = FALSE)
csvrnames = data.frame(row.names(csvproc[[1]]),stringsAsFactors = FALSE)
colnames(csvrnames) = "Parameters"

reflist = unique((sapply(csvproc,get_ref)))
masterdatalist=lapply(reflist,extract_reference,argcsvfile=csvproc,argcsvrnames=csvrnames)
masterframe = bind_rows(masterdatalist)
write.csv(masterframe,paste(sampleid,"all_coverage","csv",sep="."), na = "NO DATA", row.names = FALSE, quote = FALSE)
