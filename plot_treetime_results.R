rm(list=ls())

require(ape)
require(ggplot2)
require(ggtree)

# read in the .nexus timetree generated from the iq-tree kimura model
# as well as the data file from Ben S. with lineage info

#read in the timed tree:
#tree <- read.nexus('constant_site_correction_kimura/treetime_results/timetree.nexus')

# read in the un-timed tree:
tree <- read.tree('constant_site_correction_kimura/Malawi_final_filtered.treefile')

tree1 <- read.tree('lineage1_tree.newick')
tree2 <- read.tree('lineage2_tree.newick')
tree3 <- read.tree('lineage3_tree.newick')
tree4 <- read.tree('lineage4_tree.newick')
treeMbovis <- read.tree('lineageMbovis_tree.newick')

timetree1 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage1_results/timetree.nexus')
timetree2 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage2_results/timetree.nexus')
timetree3 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage3_results/timetree.nexus')
timetree4 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage4_results/timetree.nexus')
timetreeMbovis <- ape::read.nexus('constant_site_correction_kimura/treetime_lineageMbovis_results/timetree.nexus')


# ensure the tree is a bifurcating, rooted tree
if (!is.binary(tree)){# if not bifuricating tree
  tree=multi2di(tree) # Change to bifuricating
  tree$edge.length[which(tree$edge.length==0)]<-0.000000000001 # make all branch length non-zero
}

drugdat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')

#merge the tree and drugdat
drugdat.merge <- merge(data.frame(Sequence.name=tree$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )

drugdat1.merge <- merge(data.frame(Sequence.name=timetree1$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )

drugdat2.merge <- merge(data.frame(Sequence.name=timetree2$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )

drugdat3.merge <- merge(data.frame(Sequence.name=timetree3$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )

drugdat4.merge <- merge(data.frame(Sequence.name=timetree4$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )

drugdatMbovis.merge <- merge(data.frame(Sequence.name=timetreeMbovis$tip.label), drugdat,
	by.x="Sequence.name",by.y="Sequence.name",all.x=TRUE )


#merge the tree and blastdat
blastdat.merge <- merge(data.frame(Sequence_name=tree$tip.label), blastdat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

# Create a ggtree object for drugdat
p.drugdat <- ggtree(tree,layout='circular')

# lineage-specific ggtree objects:
p1.drugdat <- ggtree(timetree1,layout='circular')
p2.drugdat <- ggtree(timetree2,layout='circular')
p3.drugdat <- ggtree(timetree3,layout='circular')
p4.drugdat <- ggtree(timetree4,layout='circular')
pMbovis.drugdat <- ggtree(timetreeMbovis,layout='circular')



# Join the drugdat data to the ggtree object
p.drugdat <- p.drugdat %<+% drugdat.merge

p1.drugdat <- p1.drugdat %<+% drugdat1.merge
p2.drugdat <- p2.drugdat %<+% drugdat2.merge
p3.drugdat <- p3.drugdat %<+% drugdat3.merge
p4.drugdat <- p4.drugdat %<+% drugdat4.merge
pMbovis.drugdat <- pMbovis.drugdat %<+% drugdatMbovis.merge

# Color the tree based on a metadata column, for example, 'Major.Lineage'

#looks like iqtree2 recovered the lineages:
p.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()

p1.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()

p2.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()

p3.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()

p4.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()

pMbovis.drugdat + geom_tippoint(aes(color=Major.Lineage)) +
    theme_tree2()



#look at some drug resistance data:
p.drugdat + geom_tippoint(aes(color=Drug.resistance.Tbprofiler)) +
    theme_tree2()


# Create a ggtree object for blastdat
p.blastdat <- ggtree(tree,layout='circular')

# Join the blastdat data to the ggtree object
p.blastdat <- p.blastdat %<+% blastdat.merge
 
# Color the tree based on a metadata column, for example, 'Major.Lineage'
p.blastdat + geom_tippoint(aes(color=wardID)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=x20hiv)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=x17tbtype)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=sympcough)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=sympweight)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=coughduration)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=radio)) +
    theme_tree2()

p.blastdat + geom_tippoint(aes(color=phone)) +
    theme_tree2()









## this stuff below was all using the ace function, which is throwing some
# incomprehensible error on my desktop. Just drawing metadata directly onto
# the untimed tree now 

lineage <- dat$Major.Lineage
names(lineage) <- dat$Sequence.name

lineage = lineage[which(names(lineage) %in% tree$tip.label)]

#this is now throwing an error on the desktop. Not sure why.
#anc_lineage = suppressWarnings(ace(lineage,tree,type='d',method='ML'))

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


