#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
SAMPLE = args[2]
pandocloc = args[3]

library(ggplot2)
library(dplyr)
library(tinytex)

mastervector = readLines(inputfile)
library(data.table)
fields = lapply(mastervector, function(x,split){unlist(strsplit(x,split))}, split="\t")
fields = as.data.frame(data.table::transpose(fields))
names(fields)<-c("Chromosome", "Position","Coverage")
all_coverage_tbl <- fields %>%
  group_split(Chromosome)

i = 1

pltlist = list()

for (coveragetbl in all_coverage_tbl) {
  coveragetbl$Position <-as.numeric(levels(coveragetbl$Position))[coveragetbl$Position]
  coveragetbl$Coverage <-as.numeric(levels(coveragetbl$Coverage))[coveragetbl$Coverage]
  coveragetbl$Chromosome <-as.character(levels(coveragetbl$Chromosome))[coveragetbl$Chromosome]
  plt = ggplot(data=coveragetbl, aes(x=Position, y=Coverage, group=1)) + geom_line() + ggtitle(coveragetbl[1,]$Chromosome)+ ylim(0,NA)
  pltlist[[i]] = plt 
  i = i+1
}

markdownfile = paste(SAMPLE,".","basedepth.rmd",sep="")

Sys.setenv(R_PANDOC = pandocloc)
if (file.exists(markdownfile)) 
{file.remove(markdownfile)}

write(paste('---'),file = markdownfile,append=TRUE)
write(paste('header-includes:'),file = markdownfile,append=TRUE)
write(paste('- \\usepackage{multicol}'),file = markdownfile,append=TRUE)
write(paste('- \\newcommand{\\btwocol}{\\begin{multicols}{2}}'),file = markdownfile,append=TRUE)
write(paste('- \\newcommand{\\etwocol}{\\end{multicols}}'),file = markdownfile,append=TRUE)

write(paste('output: pdf_document'),file = markdownfile,append=TRUE)
write(paste('---','\n'),file = markdownfile,append=TRUE)

write(paste('**Sample Name:** ',SAMPLE),file = markdownfile,append=TRUE)

write(paste('\n\\btwocol'),file = markdownfile,append=TRUE)
for(j in c(1:length(pltlist)))
{
  write(paste('```{r echo=FALSE, fig.width=5, fig.height=3}','\n','print(pltlist[[',j,']])','\n','```',sep=''),file = markdownfile,append=TRUE)
}

write(paste('\n\\etwocol'),file = markdownfile,append=TRUE)

rmarkdown::render(markdownfile)
