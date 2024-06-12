rm(list=ls())

require(ape)
require(ggplot2)
require(ggtree)

# read in the .nexus timetree generated from the iq-tree kimura model
# as well as the data file from Ben S. with lineage info

tree <- read.nexus('timetree.nexus')

# ensure the tree is a bifurcating, rooted tree
if (!is.binary(tree)){# if not bifuricating tree
  tree=multi2di(tree) # Change to bifuricating
  tree$edge.length[which(tree$edge.length==0)]<-0.000000000001 # make all branch length non-zero
}


dat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')

lineage <- dat$Major.Lineage
names(lineage) <- dat$Sequence.name

lineage = lineage[which(names(lineage) %in% tree$tip.label)]

anc_lineage = suppressWarnings(ace(lineage,tree,type='d',method='ML'))

node_anc <- colnames(anc_lineage$lik.anc)[apply(anc_lineage$lik.anc, 1, which.max)]

names(node_anc)<-row.names(anc_lineage$lik.anc)

node_anc_colors <- node_anc
node_anc_colors[which(node_anc_colors=="lineage1")]<-"violet"
node_anc_colors[which(node_anc_colors=="lineage2")]<-"lightblue"
node_anc_colors[which(node_anc_colors=="lineage3")]<-"palegreen"
node_anc_colors[which(node_anc_colors=="lineage4")]<-"salmon1"
node_anc_colors[which(node_anc_colors=="M.bovis")]<-"gold"


lineageplot <- ggtree(tree,right=TRUE, mrsd = max(dat$date),color="lightgrey") + 
  geom_nodepoint(color=node_anc_colors) + theme_tree2()

# the major clades appear to be the main lineages! 

#look at other characteristics:

resistance <- dat$Drug.resistance.Tbprofiler
names(resistance) <- dat$Sequence.name

resistance = resistance[which(names(resistance) %in% tree$tip.label)]

anc_resistance = suppressWarnings(ace(resistance,tree,type='d',method='ML'))

node_anc <- colnames(anc_resistance$lik.anc)[apply(anc_resistance$lik.anc, 1, which.max)]

names(node_anc)<-row.names(anc_lineage$lik.anc)

node_anc_colors <- node_anc
node_anc_colors[which(node_anc_colors=="Drug-resistant")]<-"black"
node_anc_colors[which(node_anc_colors=="MDR")]<-"red"
node_anc_colors[which(node_anc_colors=="Sensitive")]<-"lightgrey"

drugplot <- ggtree(tree,right=TRUE, mrsd = max(dat$date),color="lightgrey") + 
  geom_nodepoint(color=node_anc_colors) + theme_tree2()

# color the BLAST data on the tree (blastdat)

library(RColorBrewer)

trait <- blastdat[,k]
names(trait) <- blastdat[,3] 

anc_trait = suppressWarnings(ace(trait,tree,type='d',method='ML'))

node_anc <- colnames(anc_trait$lik.anc)[apply(anc_trait$lik.anc, 1, which.max)]

names(node_anc)<-row.names(anc_trait$lik.anc)
node_anc_colors <- node_anc
cols = brewer.pal(length(unique(node_anc_colors)), "BrBG")
for(i in 1:length(unique(node_anc_colors)) ){
	node_anc_colors[which(node_anc_colors==unique(node_anc_colors)[i])]<-cols[i]
}

traitplot <- ggtree(tree,right=TRUE, mrsd = max(dat$date),color="lightgrey") + 
  geom_nodepoint(color=node_anc_colors) + theme_tree2()

# itol: interactive tree of life
# figtree - try to map on the metadata
# bisse models
# phylogenetic signal, ideally without time and without inferred node states
# try circular phylogeny
# local branching index! probably needs a time tree - within the individual clades? could try it with
#  untimed tree
# think about using PCA for dimensionality reduction, or formulate a lasso regression to predict local branching index

# 1. PCA dimensionality in the metadata
# 2. What is ancestral state reconstruction, anyway?
# 3. local branching index in the timed trees, seeing if PCA-ish outputs are useful predictors of LBI in a regression
# 4. Do any of the metadata variables appear significant in a regression predicting local branching index?


