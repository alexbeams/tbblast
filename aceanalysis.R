rm(list=ls())

# Calculate sequence densities using patristic distances weighted by decaying exponentials
#	(similar to THD and LBI)
 
library(ape)
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

#keep track of the sequences in each lineage
lin1nms <- tree1$tip.label
lin2nms <- tree2$tip.label
lin3nms <- tree3$tip.label
lin4nms <- tree4$tip.label
linMbovnms <- treeMbovis$tip.label

timetree1 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage1_results/timetree.nexus')
timetree2 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage2_results/timetree.nexus')
timetree3 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage3_results/timetree.nexus')
timetree4 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage4_results/timetree.nexus')
timetreeMbovis <- ape::read.nexus('constant_site_correction_kimura/treetime_lineageMbovis_results/timetree.nexus')

# load in the metadata and the lineage/drug resistance data
drugdat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')
names(drugdat)[1] <- 'Sequence_name'

# merge them by sequence name
dat <- merge(blastdat, drugdat, by='Sequence_name')


# ancestral character estimation

#x04fac_code

lin1inds <- which(dat$Sequence_name %in% tree1$tip.label)
x <- dat[lin1inds,'x04fac_code']
names(x) <- dat[lin1inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree1, type='d')


co <- sample(colors(), length(levels(x)))


plot(tree1, type = "tidy", FALSE, label.offset = 0,cex=0.2,main='Lineage 1: x04fac_code')
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)
legend('bottomleft',pch=19,col=co,legend=colnames(ans$lik.anc),bty='n')

lin2inds <- which(dat$Sequence_name %in% tree2$tip.label)
x <- dat[lin2inds,'x04fac_code']
names(x) <- dat[lin2inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree2, type='d')

plot(tree2, type = "tidy", FALSE, label.offset = 0,cex=0.5)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.5)

lin3inds <- which(dat$Sequence_name %in% tree3$tip.label)
x <- dat[lin3inds,'x04fac_code']
names(x) <- dat[lin3inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree3, type='d')

plot(tree3, type = "tidy", FALSE, label.offset = 0,cex=0.3)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

lin4inds <- which(dat$Sequence_name %in% tree4$tip.label)
x <- dat[lin4inds,'x04fac_code']
names(x) <- dat[lin4inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree4, type='d')

plot(tree4, type = "tidy", FALSE, label.offset = 0,cex=0.1)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.1)

# economic variables:
# let's try something with fewer states, like radio

# phone
co <- c('blue','red')

lin1inds <- which(dat$Sequence_name %in% tree1$tip.label)
x <- dat[lin1inds,'phone']
names(x) <- dat[lin1inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree1, type='d')

plot(tree1, type = "tidy", FALSE, label.offset = 0,cex=0.2,main='Lineage 1: phone')
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)
legend('bottomleft',pch=20,col=co,legend=colnames(ans$lik.anc),bty='n',title='Phone')



lin2inds <- which(dat$Sequence_name %in% tree2$tip.label)
x <- dat[lin2inds,'phone']
names(x) <- dat[lin2inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree2, type='d')

plot(tree2, type = "tidy", FALSE, label.offset = 0,cex=0.5)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.5)

lin3inds <- which(dat$Sequence_name %in% tree3$tip.label)
x <- dat[lin3inds,'phone']
names(x) <- dat[lin3inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree3, type='d')

plot(tree3, type = "tidy", FALSE, label.offset = 0,cex=0.3)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

lin4inds <- which(dat$Sequence_name %in% tree4$tip.label)
x <- dat[lin4inds,'phone']
names(x) <- dat[lin4inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree4, type='d')

plot(tree4, type = "tidy", FALSE, label.offset = 0,cex=0.1)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

# radio

lin4inds <- which(dat$Sequence_name %in% tree4$tip.label)
x <- dat[lin4inds,'radio']
names(x) <- dat[lin4inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree4, type='d')

plot(tree4, type = "tidy", FALSE, label.offset = 0,cex=0.1)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

# hivstatus variables:

# x20hiv
co <- c('blue','red')

# something about the distribution of HIV positives among lineage 1
# makes the method want to attribute a strong phylogenetic signal to
# x20hiv
lin1inds <- which(dat$Sequence_name %in% tree1$tip.label)
x <- dat[lin1inds,'x20hiv']
names(x) <- dat[lin1inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree1, type='d')

plot(tree1, type = "fan", FALSE, label.offset = 0,cex=0.1,
	main='Lineage 1: hivstatus')
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = .8)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.4)


lin2inds <- which(dat$Sequence_name %in% tree2$tip.label)
x <- dat[lin2inds,'x20hiv']
names(x) <- dat[lin2inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree2, type='d')

plot(tree2, type = "tidy", FALSE, label.offset = 0,cex=0.5)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.5)

lin3inds <- which(dat$Sequence_name %in% tree3$tip.label)
x <- dat[lin3inds,'x20hiv']
names(x) <- dat[lin3inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree3, type='d')

plot(tree3, type = "tidy", FALSE, label.offset = 0,cex=0.3)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

lin4inds <- which(dat$Sequence_name %in% tree4$tip.label)
x <- dat[lin4inds,'x20hiv']
names(x) <- dat[lin4inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree4, type='d')

plot(tree4, type = "tidy", FALSE, label.offset = 0,cex=0.1)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)


# x16patcat variables:

co <- c('blue','red', 'darkgreen', 'orange', 'purple')

lin1inds <- which(dat$Sequence_name %in% tree1$tip.label)
x <- dat[lin1inds,'x16patcat']
names(x) <- dat[lin1inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree1, type='d')

plot(tree1, type = "tidy", FALSE, label.offset = 0,cex=0.2)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.2)


lin2inds <- which(dat$Sequence_name %in% tree2$tip.label)
x <- dat[lin2inds,'x16patcat']
names(x) <- dat[lin2inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree2, type='d')

plot(tree2, type = "tidy", FALSE, label.offset = 0,cex=0.5)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.5)

lin3inds <- which(dat$Sequence_name %in% tree3$tip.label)
x <- dat[lin3inds,'x16patcat']
names(x) <- dat[lin3inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree3, type='d')

plot(tree3, type = "tidy", FALSE, label.offset = 0,cex=0.3)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)

lin4inds <- which(dat$Sequence_name %in% tree4$tip.label)
x <- dat[lin4inds,'x16patcat']
names(x) <- dat[lin4inds,'Sequence_name']
x <- factor(x)
ans <- ace(x, tree4, type='d')

plot(tree4, type = "tidy", FALSE, label.offset = 0,cex=0.1)
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.3)


