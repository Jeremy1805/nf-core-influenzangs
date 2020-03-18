#!/usr/bin/env Rscript

#SECTION I: Functions for identifying good enough references and determining pass/fail condition

get_data<-function(filename){
  splitfilename=unlist(strsplit(filename,'[.]'))
  read.csv(filename, header=TRUE, sep="\t", row.names = 1, col.names = c("buff",splitfilename[2]))
}

rank_cons <- function(sampleid,samplesegment,covext,FILTER_ENV_ARG) {
  #Get list of files with the correct extension 
  file_list<-list.files(pattern=paste(sampleid,".*_",samplesegment,".",covext,sep=""))
  
  #Collect data from the coverage summary files and save into a list of data frames
  master_csv<-lapply(file_list,get_data)
  
  #Collapse the data frames into a single data frame
  mastable = bind_cols(master_csv)
  row.names(mastable)<-row.names(master_csv[[1]])
  
  #Collect information on the coverage breadth at 10% coverage depth for samples 
  tenp_vector = as.numeric(sub("%", "",as.character(unlist(mastable["10x",]))))/100
  maximal_tenp = max(tenp_vector)
  
  mastnames = names(mastable)
  
  if (maximal_tenp >= 0.9) {
    #For a sample that passes (at least one reference at >90% coverage) pick only those samples with decent 10% coverage
    chosenone = mastnames[which(tenp_vector>=0.80)]
    FILTER_ENV_ARG$state = "PASS"
  } else {
    FILTER_ENV_ARG$state = "FAIL"
    chosenone=mastnames
  }
  
  return(chosenone)
}

#SECTION II: Functions for determining type and subtype mixes in samples

cov_get_similarity <- function(test_ref,hypothesis_ref,sampleid_arg,covext_arg) {
  if (test_ref == hypothesis_ref) {
    return(0)
  } else {
  #PARAMETERS
    hypothesis_maximum = 0
    min_reads = 10
  
    test_cov = paste(sampleid_arg,hypothesis_ref,covext_arg,sep=".")
    hypothesis_cov = paste(test_ref,hypothesis_ref,covext_arg,sep=".")

    coverage_table_test=read.csv(test_cov,header=FALSE,sep="\t")
    names(coverage_table_test) = c("CHROM","POS","READS")
    coverage_table_hypothesis = read.csv(hypothesis_cov,header=FALSE,sep="\t")
    names(coverage_table_hypothesis) = c("CHROM","POS","READS")

    Relevant_pos = coverage_table_hypothesis[which(coverage_table_hypothesis$READS <= hypothesis_maximum),]$POS 
    no_intersect_positions= which((coverage_table_test$POS %in% Relevant_pos) & (coverage_table_test$READS >= min_reads))

    Percentage_positions=length(no_intersect_positions)/nrow(coverage_table_hypothesis)

    return(Percentage_positions)
  }
}

single_row_get_similarity <- function(hypothesis_arg,testvector,id_arg,covext_a) {
  return(sapply(testvector,cov_get_similarity,hypothesis_ref=hypothesis_arg,sampleid_arg=id_arg,covext_arg=covext_a))
}

powset <- function (set) { 
  n <- length(set)
  masks <- 2^(1:n-1)
  sapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] )
}

getmatch <- function (querymem,equivalentls) {
  sum(as.integer(sapply(equivalentls,FUN=match,x=querymem,nomatch=0) > 0))
}

matchall <- function (equivalencelist,membervc)  {
  return (sum(abs(sapply(membervc,FUN=getmatch,equivalentls=equivalencelist)-1)))
}

getiniteq <- function (colin,membervec) {
  membervec[which(colin)]
}

getequivalencesets <- function(filtered_references,sampleid,covext,FILTER_ENV_ARG){
  #EACH COLUMN HAS THE SAME TEST REFERNCE, EACH ROW HAS THE SAME HYPOTHESIS REFERENCE
  eval_matrix = sapply(filtered_references,single_row_get_similarity, testvector=filtered_references, id_arg = sampleid, covext_a=covext)
  testmatrix = eval_matrix<0.20
  boolmatrix = testmatrix|t(testmatrix)
  
  #get all unique equivalence classes. 
  #This is a list, where each entry in the list is a vector of references within the equivalence class. 
  list_of_eq = alply(boolmatrix,2,getiniteq,membervec=filtered_references)
  list_of_eq=unique(list_of_eq)
  
  #Create a power set of equivalence classes
  powerlist=powset(list_of_eq)
  powerlist[[1]]=NULL
  
  #Resolve any ambiguities. Penalty_vec assigns to each possible equivalence class a 
  penalty_vec=sapply(powerlist,matchall,membervc=filtered_references)
  min_penalty_list = powerlist[[match(min(penalty_vec),penalty_vec)]]
  
  if(length(min_penalty_list)>1) {
    FILTER_ENV_ARG$mix = "MIX"
  } else {
    FILTER_ENV_ARG$mix = "PASS"
  }
  
  return(min_penalty_list)
}

#SECTION III: Functions for re-ranking within each mix 

rerank_cons <- function (ref_list, sampleid_a, covsumext_a) {
  file_list <- paste(sampleid_a,ref_list,covsumext_a,sep=".")
  master_csv<-lapply(file_list,get_data)
  
  mastable = bind_cols(master_csv)
  row.names(mastable)<-row.names(master_csv[[1]])
  
  genome_size = as.numeric(as.character(unlist(mastable["Total Genome Size",])))
  zero_vector = as.numeric(as.character(unlist(mastable["Zero Coverage",])))
  nonzero_percent = (1-zero_vector/genome_size)*100
  
  tenx_vector = as.numeric(sub("%", "",as.character(unlist(mastable["10x",]))))
  hundredx_vector = as.numeric(sub("%", "",as.character(unlist(mastable["100x",]))))
  
  weighted_sum = 0.25*nonzero_percent + 0.35*tenx_vector + 0.4*hundredx_vector
  #chosenone=names(mastable[,which(zero_<=minimal_zero+0.1)])
  
  orderedsum = order(-weighted_sum)
  
  if (length(orderedsum)<3) {
    finalsize = length(orderedsum)
  } else {
    finalsize = 3
  }
  
  mastnames = names(mastable)
  
  chosenone = sort(mastnames[orderedsum[1:finalsize]])
  return(chosenone)
}

#SECTION IV: Functions for choosing the best reference 

add_col <- function(vec,namestr) {
  mat = matrix(nrow=length(vec),ncol=2)
  mat[,1] = namestr
  mat[,2] = vec
  return(mat)
} 

write_table_ls <-function(namedvec,inputnm) {
  write.table(namedvec, file = paste(inputnm,namedvec[1,1],"posmatrix",sep="."),row.names=FALSE, col.names=FALSE, na="", sep="\t", quote=FALSE)
}

create_reference_tables <- function (inputfile){
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
}

get_percentage_similarity <- function(refseg, id_arg, tupleid_arg) {
  return(as.numeric(system2("get_perc_ref.sh",c(id_arg,refseg, tupleid_arg),stdout=TRUE)))
}

get_best_reference <- function (chosenones, sampleid_arg){
  if(length(chosenones) > 1) {
    top_reflist=paste(chosenones,collapse=" ")
    tupleid=paste(chosenones,collapse="")
  
    clustal_align_command=paste("\"",top_reflist,"\" ",tupleid,sep="")
    system2("clustal_align.sh",args=c(paste("\"",top_reflist,"\"",sep=""),tupleid))
  
    create_reference_tables(paste(tupleid,"fasta",sep="."))
    percentage_coverage = sapply(chosenones, get_percentage_similarity, id_arg=sampleid, tupleid_arg=tupleid)
  
    return(chosenones[which.max(percentage_coverage)])
    
  } else {
  return(chosenones)  
  }
}

library("Biostrings")
library("dplyr")
library("plyr")

args = commandArgs(trailingOnly=TRUE)
sampleid = args[1]
samplesegment = args[2]

covsumext="cov.summary"
covext="cov"

FILTER_STATE_ENV <- new.env()

filtered_references = rank_cons(sampleid,samplesegment,covsumext,FILTER_STATE_ENV)
init_filter = FILTER_STATE_ENV$state

if (length(filtered_references) == 1) {
  best_ref = filtered_references
} else if(FILTER_STATE_ENV$state == "PASS") {
  equivalent_references = getequivalencesets(filtered_references,sampleid,covext,FILTER_STATE_ENV)

  top_3_references = lapply(equivalent_references, rerank_cons, sampleid_a=sampleid, covsumext_a=covsumext)

  best_ref = sapply(top_3_references, get_best_reference, sampleid_arg=sampleid)
} else {
  
  top_3_references = lapply(list(filtered_references), rerank_cons, sampleid_a=sampleid, covsumext_a=covsumext)
  
  best_ref = sapply(top_3_references, get_best_reference, sampleid_arg=sampleid)
}

if ( (length(FILTER_STATE_ENV$mix) == 0) | FILTER_STATE_ENV$state == "FAIL"){
  filter_output = paste(FILTER_STATE_ENV$state,paste(best_ref,collapse="/"),sep=":")
} else {
  filter_output = paste(FILTER_STATE_ENV$mix,paste(best_ref,collapse="/"),sep=":")
}

best_ref_output=paste(paste(sampleid,best_ref,init_filter,sep="."),collapse="\n")
write(filter_output,paste(sampleid,samplesegment,"filter",sep="."))
write(best_ref_output,paste(sampleid,samplesegment,"finalchosen",sep="."))