#!/usr/bin/env Rscript

get_id_vec<-function(filename){
    splitfilename=unlist(strsplit(filename,'[.]'))
    return(splitfilename[1])
}

get_seg_vec<-function(filename){
  splitfilename=unlist(strsplit(filename,'[.]'))
  return(splitfilename[2])
}

get_state<-function(seg,id){
  #Construct the expected filenames from id and segment
  filterfile=paste(id,seg,"filter",sep=".")
  #Initial filter:"Pass" if more than 90% coverage in one reference or "Fail" otherwise
  #If "Pass" the first filter, then check for mixed samples:either "Mix" or "Pass"
  #Filterfile also includes list of best references
  if (file.exists(filterfile)) {
    filestate = readLines(filterfile) #Get mixture filter state (Pass/Fail/Mix) and list of references
    } else {
    filestate = NA    
  }
  
  return(filestate)
}

get_state_vec <- function(segvec,sampleid){
  return(sapply(segvec,get_state, id=sampleid))
}

get_state_matrix <- function(idvec,segvector){
  return(sapply(idvec,get_state_vec, segvec=segvector))
}

#Get list of files that have the filter extension. 
#This is done only to extract a list of unique IDs and unique segments  
file_list<-list.files(pattern=paste("*","filter",sep="."))

#Get segments and sample ids in vector form from the list of files
idvec = sapply(file_list,get_id_vec)
segvec = sapply(file_list,get_seg_vec)

#Extract only unique id and unique segments
uniqid = unique(idvec)
uniqseg = unique(segvec)
rowz = length(uniqid)
colz = length(uniqseg)

#Make a csv file containing pass/fail condition of filters and the chosen references.
datamatrix=t(get_state_matrix(uniqid,uniqseg))
write.csv(datamatrix,"all_samples_filtered.csv")