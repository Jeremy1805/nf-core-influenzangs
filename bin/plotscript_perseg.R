#!/usr/bin/env Rscript

get_field<- function(filename) {
  if (!file.exists(filename)){return(-1)}
  mastervector = readLines(filename)
  fields = unlist(strsplit(mastervector,"\t")) 
  if (length(fields) != 29){return(-1)}
  if (fields[2] != "Total Genome Size"){return(-1)}
  if (fields[4] != "Zero Coverage"){return(-1)}
  if (fields[6] != "Non Zero Coverage"){return(-1)}
  if (fields[8] != "Low Coverage <= 5"){return(-1)}
  if (fields[10] != "Min Coverage"){return(-1)}
  if (fields[12] != "Max Coverage"){return(-1)}
  if (fields[14] != "Median Coverage"){return(-1)}
  if (fields[16] != "Mean Coverage"){return(-1)}
  if (fields[18] != "Depth Coverage Level"){return(-1)}
  if (fields[19] != "Percentage of Kit Covered"){return(-1)}
  if (fields[20] != "10x"){return(-1)}
  if (fields[22] != "20x"){return(-1)}
  if (fields[24] != "50x"){return(-1)}
  if (fields[26] != "100x"){return(-1)}
  if (fields[28] != "1000x"){return(-1)}
  
  fnamevec = unlist(strsplit(fields[1],"[.]"))
  name = fnamevec[1]
  ref = fnamevec[2]
  
  numvec = c()
  numvec[1] = as.numeric(fields[3])
  numvec[2] = as.numeric(sub("%", "",fields[21]))

  output <- data.frame("Filename" = filename,"Sample_name" = name,
                       "Reference" = ref,"Total_Genome_Size" = numvec[1],"tenx_coverage" = numvec[2],stringsAsFactors = FALSE)
  return(output)
}

get_all_field_no_list<- function(mastervec) {
  outputvec = c()
  errorvec = c()
  masterveclen = length(mastervec)
  for(i in c(1:masterveclen))
  {
    outtemp = get_field(mastervec[i])
    if (outtemp[1] != -1)
    {outputvec = rbind(outputvec,outtemp)}
    else
    {errorvec <- rbind(errorvec,mastervec[i])}
  }
  rownames(outputvec)<-c()
  return(list(outputvec,errorvec))
}

get_ref<- function(id,seg) {
  reffile=paste(id,seg,"finalchosen",sep=".")
  chosenref = sapply(strsplit(readLines(reffile),"[.]"),"[[",2)
  return(chosenref)
}

mark_chosen_one<- function(data,segment) {
  samples = unique(data[,2]) #Gets a vector of only the id
  names(samples)=samples #names the vector
  chosenrefvec<-lapply(samples,get_ref,seg=segment)
  datalen = nrow(data)
  for(i in c(1:datalen))
  {
    nametemp = data[i,"Sample_name"]
    if(data[i,"Reference"] %in% chosenrefvec[[nametemp]])
    {
      data[i,"chosenref"] = "Chosen"
    }
    else
    {
      data[i,"chosenref"] = "Not chosen"
    }
  }
  return(data)
}

#df = read.table(args[1], header=TRUE)
library(ggplot2)
library(plotly)

args = commandArgs(trailingOnly=TRUE)
flusegment = args[1]
covsumext = args[2]

file_list<-list.files(pattern = paste("*",covsumext,sep="."))
temp = get_all_field_no_list(file_list)
all_data = temp[[1]]
all_data = mark_chosen_one(all_data,flusegment)
p=ggplot(all_data, aes(Sample_name, tenx_coverage, color=Reference, shape=factor(chosenref))) +geom_point(aes(text=sprintf("Sample name: %s<br>10x depth percent: %s", Sample_name, tenx_coverage)))+theme(axis.text.x = element_text(angle=90),legend.position="top")+labs(x="",y="Percent of Bases with \u2265 10% Depth", color = "References", shape = "")
picname = paste("10x_coverage",flusegment,sep="_")
ggsave(paste(picname,"png",sep="."),plot=p,device="png",width=40,height = 7)
pl <- ggplotly(p,tooltip=c("text","Reference"))
htmlwidgets::saveWidget(pl, paste(picname,"html",sep="."))
