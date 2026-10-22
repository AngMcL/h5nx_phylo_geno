#!/usr/bin/env Rscript
# Root IQtree output by rtt

#usage:
# Rscript root-iq-rtt.R "folder_with_trees"
library(tidyverse)
library(gtools)
library(ape)

# Input trees from fasttree
args<-commandArgs(trailingOnly=TRUE)
f.in.tre<-args[1] #folder with trees

tre.in<-list.files(f.in.tre, pattern = ".treefile$",full.names = T)  %>% mixedsort()
l.tre<-length(tre.in)

# Output
today<-Sys.Date()
f.out<-paste0(f.in.tre,"_rooted/")
if(!dir.exists(f.out)){dir.create(f.out)}

#Rooted trees out
tre.out<-str_replace_all(tre.in, f.in.tre, f.out) %>%
  str_replace_all(".fasta","iqtree_rooted")

#fill in dates if needed with sprintf for leading zero
mm<-sprintf("%02d", 1:12)
dd<-sprintf("%02d", 1:28) #up to 28 bc feb
seed<-1 #repeatable random sampling of incomplete dates

for(i in 1:l.tre){
  tree<-ape::read.tree(tre.in[i]) 
  # is.rooted(tree)
  # is.binary(tree)
  
  dates.t<-sapply(tree$tip.label, function(x) unlist(str_split(x, "\\|") )[2])
  #instances of incomplete date or month
  dates.yr<-which(sapply(dates.t, function(x) length(unlist(str_split(x, "-"))))==1) #year
  dates.yrm<-which(sapply(dates.t, function(x) length(unlist(str_split(x, "-"))))==2) #year-mm
  #sample random and concat, replace
  dates.t[dates.yr]<-sapply(dates.t[dates.yr], function(x) paste(x, sample(mm,1), sample(dd,1),sep="-"))
  dates.t[dates.yrm]<-sapply(dates.t[dates.yrm], function(x) paste(x, sample(dd,1),sep="-"))
  #make decimal date
  dates.dec<-decimal_date(as.Date(dates.t))
  tree.root<-ape::rtt(t = tree, tip.dates = dates.dec, objective = "rms")
  #check
  is.rooted(tree.root)
  #export
  write.tree(tree.root, tre.out[i])
}
