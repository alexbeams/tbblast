rm(list=ls())

#load in the time trees and color in the leaves according to the 
# characteristics in the Metadata file

require("ape")
require("phytools")
require("ggalluvial")
require(ggplot2)
require(ggtree)
require(tidyr)
require(dplyr)
require(tibble)


tree <- read.tree("timetree.nexus")

Metadata <- read.csv('~/Documents/TB/BLAST_ePAL_SeqID_NoGPS.csv')

tree <- multi2di(tree)
tree <- as.phylo(tree)

anc_Lineage <- suppressWarnings(ace(MaybeLineage,tree,type='d',method='ML'))

node_anc <- colnames(anc_Lineage$lik.anc)[apply(anc_Lineage$lik.anc, 1, which.max)]
names(node_anc)<-row.names(anc_Lineage$lik.anc)

node_anc_colors <- node_anc
node_anc_colors[which(node_anc_colors=="1")]<-"violet"
node_anc_colors[which(node_anc_colors=="2")]<-"lightblue"
node_anc_colors[which(node_anc_colors=="3")]<-"palegreen"
node_anc_colors[which(node_anc_colors=="4")]<-"salmon1"
node_anc_colors[which(node_anc_colors=="reseq")]<-"gold"


