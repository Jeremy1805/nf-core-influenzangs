#!/usr/bin/env Rscript

add_col <- function(vec,namestr) {
  mat = matrix(nrow=length(vec),ncol=2)
  mat[,1] = namestr
  mat[,2] = vec
  return(mat)
} 

write_table_ls <-function(namedvec,inputnm) {
  write.table(namedvec, file = paste(inputnm,namedvec[1,1],"posmatrix",sep="."),row.names=FALSE, col.names=FALSE, na="", sep="\t", quote=FALSE)
}

library("Biostrings")

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
inputname = unlist(strsplit(inputfile,"[.]"))[1]

masterfastaset = readDNAStringSet(inputfile)

#Check that all sizes are the same 
if(all(width(masterfastaset) == width(masterfastaset[1])) != TRUE) {
  print("ERROR: ALIGNMENT DID NOT YIELD EQUALLY SIZED SEQUENCES")
  quit(status=10)
}

samplenumber = length(masterfastaset)
baselength = width(masterfastaset[1])

#Create base location vectors
basenumvec = matrix(data=0, nrow = samplenumber, ncol = baselength)
fastacharvec = matrix(data=0, nrow = samplenumber, ncol = baselength)

#create a character matrix of bases at each alignment position for each sample
for (samplei in c(1:samplenumber)) {
  fastacharvec[samplei,] = unlist(strsplit(as.character(masterfastaset[[samplei]]),""))
}

#create an integer matrix of individual fasta positions at each alignment position for each sample
for (samplei in c(1:samplenumber)) {
  basenumvec[samplei,] = cumsum(as.integer(fastacharvec[samplei,]!="-"))
}

#Determine locations where a discrepancy exists
booleancheck = consensusMatrix(masterfastaset)==samplenumber
differencevector = !apply(booleancheck,2,any)

#Record indices where differences occur
differenceindexes = which(differencevector)

#Get bases and actual base positions at difference points
basenumvecdiff = t(basenumvec[,differenceindexes])

poslist= split(basenumvecdiff, rep(1:ncol(basenumvecdiff), each = nrow(basenumvecdiff)))
chromposlist=mapply(FUN=add_col, vec=poslist,namestr= names(masterfastaset), SIMPLIFY = FALSE)

lapply(chromposlist,write_table_ls,inputnm=inputname)
