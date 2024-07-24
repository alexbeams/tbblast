rm(list=ls())

# to do:
# generate LBI for timed and un-timed trees separately
# for timed trees, probably calculate LBI on different clades separately
# look into other phylogenetic signal metrics
# can we use BiSSE for traits like location or HIV status?
# look into timed haplotype density calculation
# BiSSE models might be ok for the 707 sequences; might be approaching the limit of feasibility for them

# THD is meant to measure how "branchy" a sequence's ancestors were
#	it basically calculates a kernel density estimate for sequences, 
#	and if testing is random, successful clades transmit many new
#	descendants in a short time, producing a higher density
#	in that region of sequence space
#	it requires a Hamming distance matrix; probably would be better to use on
#	genes, rather than nucleotides...

# LBI is accomplishing something similar. The Neher et al paper showed that it 
#	correlates reasonably well with a Brownian-Motion model for fitness
#	along lineages


# Is there a more general way to identify which parts of the tree are "hottest"?

# Cedric's idea: try using Multi-Locus Sequence Typing on the raw reads to see whether
#	the maximum likelihood tree based on SNPs is consistent with strain type


library(ape)
library(ggplot2)
library(ggtree)
library(Biostrings)


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

# itol: interactive tree of life
# figtree - try to map on the metadata

# how are metadata features distributed in the tree? can we find a phylogenetic signal?
# bisse models
# phylogenetic signal, ideally without time and without inferred node states


# are some clades more successful than others? How can we quantify this?
# local branching index! probably needs a time tree - within the individual clades? could try it with
#  untimed tree

# other node/tip statistics we could use? 
# once we have some features which we think describe clade success, we can
# think about using PCA for dimensionality reduction, or formulate a lasso regression to predict
# the local branching index, for example

# 1. PCA dimensionality in the metadata
# 2. local branching index in the timed trees, seeing if PCA-ish outputs are useful predictors of LBI in a regression
# 3. Do any of the metadata variables appear significant in a regression predicting local branching index?


#########################
# Local branching index #
#########################

# https://github.com/bdearlove/treeImbalance/
# https://rdrr.io/github/bdearlove/treeImbalance/src/R/lbi.R

#' LBI
#' 
#' This calculates the local branching index (LBI), defined by Neher, Russell and Shraiman (2014) as the the tree length surrounding a given node or tip exponentially discounted with increasing distance from that node. 
#' @param tree a tree of class \code{phylo}.
#' @param tau a number giving the scaling factor to use, suggested by Neher et al. to be 0.0625 times the average pairwise distance in the sample.
#' @param transform a function by which to transform the LBI. Default is no transformation.
#' @return an vector of class \code{numeric} giving the LBI for the tips and internal nodes in the same order as given in the tree \code{phylo} object.
#' @keywords local branching index, LBI
#' @references Neher RA, Russell CA, Shraiman BI (2014) Predicting evolution from the shape of genealogical trees. eLife 3:e03568. doi:10.7554/eLife.03568
#' @export
#' @examples
#' tree<-rtree(50)
#' tree.lbi<-lbi(tree)
#' tree$node.label<-tree.lbi[(tree$Nnode+2):length(tree.lbi)]
#' plot(tree, show.node.label = TRUE)

lbi<-function(tree, tau=0.0005,transform=function(x){x}){
  nnodes<-tree$Nnode
  ntips<-tree$Nnode+1

  node.times<-dist.nodes(tree)[(ntips+1),1:(2*ntips-1)]
    
  node.postorder<-as.numeric(names(sort(node.times,decreasing=T)))
  node.preorder<-as.numeric(names(sort(node.times)))
  
  #Initialise the arrays
  uppolarizer<-array(0,(2*ntips)) # 2*ntips - 1 
  downpolarizer<-array(0,(2*ntips)) # 2*ntips - 1 
  lbi<-array(0,(2*ntips))
  
  #Traverse tree postorder to calculate message to parents
  for(i in node.postorder[1:1020]){
    node.uppolarizer<-0
    
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    node.uppolarizer<-node.uppolarizer+sum(uppolarizer[node.children])

    node.edge<-which(tree$edge[,2]==i)
    if(length(node.edge>0)){
      bl<-tree$edge.length[node.edge]/tau
      node.uppolarizer<-node.uppolarizer*exp(-bl)
      node.uppolarizer<-node.uppolarizer+tau*(1-exp(-bl))
    }
    uppolarizer[i]<-node.uppolarizer
  }
  
  #Traverse tree preorder to calculate message to children
  for(i in node.preorder){
    node.downpolarizer<-downpolarizer[i]
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    
    for(j in node.children){
      downpolarizer[j]<-node.downpolarizer+sum(uppolarizer[setdiff(node.children,j)])
      
      child.edge<-which(tree$edge[,2]==j)
      bl<-tree$edge.length[child.edge]/tau
      
      downpolarizer[j]<-downpolarizer[j]*exp(-bl)+tau*(1-exp(-bl))
    }
  }
  
  #Calculate lbi
  for(i in node.postorder){
    node.lbi<-downpolarizer[i]
    
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    node.lbi<-node.lbi+sum(uppolarizer[node.children])
    
    lbi[i]<-transform(node.lbi)
  }
  return(lbi)
}

#generate local branching indices for the untimed trees:
x1 <- lbi(tree1,tau=.00001)
x1 <- x1[1:Ntip(tree1)]
x1nms <- tree1$tip.label
plot(x1)


x2 <- lbi(tree2)
x2 <- x2[1:Ntip(tree2)]
x2nms <- tree2$tip.label
plot(x2)

x3 <- lbi(tree3)
x3 <- x3[1:Ntip(tree3)]
x3nms <- tree3$tip.label
plot(x3)


x4 <- lbi(tree4,tau=.00001)
x4 <- x4[1:Ntip(tree4)]
x4nms <- tree4$tip.label
plot(x4)

xMbovis <- lbi(treeMbovis)
xMbovis <- xMbovis[1:Ntip(treeMbovis)]
xMbovisnms <- treeMbovis$tip.label
plot(xMbovis)

lbis <- c(x1,x2,x3,x4,xMbovis)
lbinms <- c(x1nms,x2nms,x3nms,x4nms,xMbovisnms)
lbidat <- data.frame(cbind(lbinms,lbis))
colnames(lbidat) <- c('Sequence_name','lbi')

dat = merge(blastdat,lbidat)

dat$lbi <- as.numeric(dat$lbi)

#merge the trees and dat
dat1.merge <- merge(data.frame(Sequence_name=tree1$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat2.merge <- merge(data.frame(Sequence_name=tree2$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat3.merge <- merge(data.frame(Sequence_name=tree3$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat4.merge <- merge(data.frame(Sequence_name=tree4$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

datMbovis.merge <- merge(data.frame(Sequence_name=treeMbovis$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

# lineage-specific ggtree objects:
p1.dat <- ggtree(tree1,layout='circular')
p2.dat <- ggtree(tree2,layout='circular')
p3.dat <- ggtree(tree3,layout='circular')
p4.dat <- ggtree(tree4,layout='circular')
pMbovis.dat <- ggtree(treeMbovis,layout='circular')



p1.dat <- p1.dat %<+% dat1.merge
p2.dat <- p2.dat %<+% dat2.merge
p3.dat <- p3.dat %<+% dat3.merge
p4.dat <- p4.dat %<+% dat4.merge
pMbovis.dat <- pMbovis.dat %<+% datMbovis.merge

# Color the tree based on a metadata column, for example, 'lbi'

p1lbi <- p1.dat + geom_tippoint(aes(color=lbi)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2() + ggtitle('lineage 1')

p2lbi <- p2.dat + geom_tippoint(aes(color=lbi)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 2')

p3lbi <- p3.dat + geom_tippoint(aes(color=lbi)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 3')

p4lbi <- p4.dat + geom_tippoint(aes(color=lbi)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 4')

pMbovislbi <- pMbovis.dat + geom_tippoint(aes(color=lbi)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage M.bovis')


#generate local branching indices for the timed trees:
xt1 <- lbi(timetree1,tau=.00001)
xt1 <- xt1[1:Ntip(timetree1)]
xt1nms <- timetree1$tip.label
plot(xt1)


xt2 <- lbi(timetree2)
xt2 <- xt2[1:Ntip(timetree2)]
xt2nms <- timetree2$tip.label
plot(xt2)

xt3 <- lbi(timetree3)
xt3 <- xt3[1:Ntip(timetree3)]
xt3nms <- timetree3$tip.label
plot(xt3)


xt4 <- lbi(timetree4,tau=.00001)
xt4 <- xt4[1:Ntip(timetree4)]
xt4nms <- timetree4$tip.label
plot(xt4)

xtMbovis <- lbi(timetreeMbovis)
xtMbovis <- xtMbovis[1:Ntip(timetreeMbovis)]
xtMbovisnms <- timetreeMbovis$tip.label
plot(xtMbovis)

lbits <- c(xt1,xt2,xt3,xt4,xtMbovis)
lbitnms <- c(xt1nms,xt2nms,xt3nms,xt4nms,xtMbovisnms)
lbitdat <- data.frame(cbind(lbinms,lbits))
colnames(lbitdat) <- c('Sequence_name','lbit')


dat = merge(dat,lbitdat)

dat$lbit <- as.numeric(dat$lbit)

#merge the trees and dat 
timedat1.merge <- merge(data.frame(Sequence_name=timetree1$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

timedat2.merge <- merge(data.frame(Sequence_name=timetree2$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

timedat3.merge <- merge(data.frame(Sequence_name=timetree3$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

timedat4.merge <- merge(data.frame(Sequence_name=timetree4$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

timedatMbovis.merge <- merge(data.frame(Sequence_name=timetreeMbovis$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

# lineage-specific ggtree objects:
p1.timedat <- ggtree(timetree1,layout='circular')
p2.timedat <- ggtree(timetree2,layout='circular')
p3.timedat <- ggtree(timetree3,layout='circular')
p4.timedat <- ggtree(timetree4,layout='circular')
pMbovis.timedat <- ggtree(timetreeMbovis,layout='circular')

p1.timedat <- p1.timedat %<+% timedat1.merge
p2.timedat <- p2.timedat %<+% timedat2.merge
p3.timedat <- p3.timedat %<+% timedat3.merge
p4.timedat <- p4.timedat %<+% timedat4.merge
pMbovis.timedat <- pMbovis.timedat %<+% timedatMbovis.merge

# Color the tree based on a metadata column, for example, 'lbit'

p1lbit <- p1.timedat + geom_tippoint(aes(color=lbit)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2() + ggtitle('lineage 1')

p2lbit <- p2.timedat + geom_tippoint(aes(color=lbit)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 2')

p3lbit <- p3.timedat + geom_tippoint(aes(color=lbit)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 3')

p4lbit <- p4.timedat + geom_tippoint(aes(color=lbit)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage 4')

pMbovislbit <- pMbovis.timedat + geom_tippoint(aes(color=lbit)) + 
	scale_color_gradient(high='red',low='blue') +
    theme_tree2()+ ggtitle('lineage M.bovis')



# Calculate the Timed Haplotype Density for the trees
#	start with tree4 and timetree4 for now


#load in the functions to calculate the THD
source('thd.R')

#create a hamming distance matrix and plug into the KDE
fasta_file <- "Malawi_final_filtered.fasta"
alignment <- readDNAStringSet(fasta_file)

compute_hamming_distance <- function(seq1, seq2) {
  seq1_chars <- as.character(unlist(strsplit(as.character(seq1), "")))
  seq2_chars <- as.character(unlist(strsplit(as.character(seq2), "")))
  sum(seq1_chars != seq2_chars)
}

# Get the number of sequences
n <- length(alignment)

# Initialize the distance matrix
hamming_matrix <- matrix(0, n, n)

# Uncomment to fill the matrix
#for (i in 1:n) {
#  for (j in i:n) {
#    dist <- compute_hamming_distance(alignment[[i]], alignment[[j]])
#    hamming_matrix[i, j] <- dist
#    hamming_matrix[j, i] <- dist
#  }
#}

load("hamming_matrix.Rdata")

m <- length(alignment[[1]])
mu <- 5e-4
timescale <- 2 
bandwidth <- tmrca2bandwidth(timescale, m , mu)
THD <- thd(hamming_matrix, bandwidth, m)

# Display hclust tree and THDs
hc <- hclust(as.dist(hamming_matrix))

par(mfrow = c(2,1))
par(mar = c(1, 5, 1, 5))
plot(hc, labels = FALSE, xlab = "", main = "", sub = "")
barplot(THD[hc$order], ylim = rev(range(THD)), ylab = "THD")


#plot with color
colors <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THD, breaks = 100))]
ph <- as.phylo(hc)
par(mfrow=c(1,1))
plot(ph,type='fan',show.tip.label=F,main='THD for 707 Sequences')
tiplabels(pch=19,col=colors,cex=1)

#probably want to run THD for each lineage separately. Just looks like
#	it's picking up lineage 4. Or, maybe just tuning the bandwidth
#	will do the job...

rownames(hamming_matrix) <- names(alignment)
colnames(hamming_matrix) <- names(alignment)

lin1inds <- sapply(lin1nms, function(x) which(names(alignment) == x  ))
lin2inds <- sapply(lin2nms, function(x) which(names(alignment) == x  ))
lin3inds <- sapply(lin3nms, function(x) which(names(alignment) == x  ))
lin4inds <- sapply(lin4nms, function(x) which(names(alignment) == x  ))
linMbovinds <- sapply(linMbovnms, function(x) which(names(alignment) == x  ))

hamming_lin1 <- hamming_matrix[lin1inds,lin1inds]
hamming_lin2 <- hamming_matrix[lin2inds,lin2inds]
hamming_lin3 <- hamming_matrix[lin3inds,lin3inds]
hamming_lin4 <- hamming_matrix[lin4inds,lin4inds]
hamming_linMbov <- hamming_matrix[linMbovinds,linMbovinds]

mu1 <- 5e-4
timescale1 <- 2 
bandwidth1 <- tmrca2bandwidth(timescale1, m , mu1)

mu2 <- 5e-4
timescale2 <- 2 
bandwidth2 <- tmrca2bandwidth(timescale2, m , mu2)

mu3 <- 5e-4
timescale3 <- 2 
bandwidth3 <- tmrca2bandwidth(timescale3, m , mu3)

mu4 <- 5e-4
timescale4 <-2 
bandwidth4 <- tmrca2bandwidth(timescale4, m , mu4)

THD1 <- thd(hamming_lin1, bandwidth1, m)
THD2 <- thd(hamming_lin2, bandwidth2, m)
THD3 <- thd(hamming_lin3, bandwidth3, m)
THD4 <- thd(hamming_lin4, bandwidth4, m)
THDMbov <- thd(hamming_linMbov, bandwidth, m)

# plot the THDs:
par(mfrow=c(2,2))

plot(THD1)
plot(THD2)
plot(THD3)
plot(THD4)




# Display hclust trees and THDs
hc1 <- hclust(as.dist(hamming_lin1))
hc2 <- hclust(as.dist(hamming_lin2))
hc3 <- hclust(as.dist(hamming_lin3))
hc4 <- hclust(as.dist(hamming_lin4))
hcMbov <- hclust(as.dist(hamming_linMbov))

#plot with color
colors1 <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THD1, breaks = 100))]
colors2 <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THD2, breaks = 100))]
colors3 <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THD3, breaks = 100))]
colors4 <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THD4, breaks = 100))]
colorsMbov <- colorRampPalette(c("blue", "red"))(100)[as.numeric(cut(THDMbov, breaks = 100))]

ph1 <- as.phylo(hc1)
ph2 <- as.phylo(hc2)
ph3 <- as.phylo(hc3)
ph4 <- as.phylo(hc4)
phMbov <- as.phylo(hcMbov)

par(mfrow=c(2,2))
plot(ph1,type='fan',show.tip.label=F,main='THD for Lineage 1 Sequences')
tiplabels(pch=19,col=colors1,cex=1)

plot(ph2,type='fan',show.tip.label=F,main='THD for Lineage 2 Sequences')
tiplabels(pch=19,col=colors2,cex=1)

plot(ph3,type='fan',show.tip.label=F,main='THD for Lineage 3 Sequences')
tiplabels(pch=19,col=colors3,cex=1)

plot(ph4,type='fan',show.tip.label=F,main='THD for Lineage 4 Sequences')
tiplabels(pch=19,col=colors4,cex=1)



# examine the ClusterPicker tree:

# Read the tree
tree <- read.tree("/Users/abeams/Documents/TB/constant_site_correction_kimura_bootstrapped/clusterpickertree.newick")

# Example: Extracting cluster names from tip labels
# Assuming tip labels are in the format "tipname_clustername"
tip_labels <- tree$tip.label
clusters <- sapply(strsplit(tip_labels, "_"), function(x) x[1])

# Create a data frame for ggtree
df <- data.frame(label = tip_labels, cluster = clusters)

# Generate a distinct color palette for 62 clusters
num_clusters <- length(unique(clusters))
colors <- brewer.pal(min(num_clusters, 12), "Set3") # Use "Set3" palette for up to 12 clusters

# If more than 12 clusters, use colorRampPalette to extend the palette
if (num_clusters > 12) {
  colors <- colorRampPalette(colors)(num_clusters)
}

# Create the ggtree plot with a fan layout
p <- ggtree(tree, layout = "fan") %<+% df + 
  geom_tippoint(aes(color = cluster), size = 2) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = setNames(colors, unique(clusters)))

# Display the plot
print(p)








# are any metadata variables correlated with LBI?

#any noticeable 1trends in time?
plot(lbi~as.factor(x05dd),dat)
plot(lbi~as.factor(x05month),dat)
plot(lbi~as.factor(x05year),dat)

#location?
plot(lbi~as.factor(wardID),dat)
plot(lbi~as.factor(x04fac_code),dat)
# what's going at the Limbe Central Ward? (is that the one?)

# age:
plot(lbi~as.factor(x12age),dat)
plot(lbi~as.factor(x13agegp),dat)

#cgh?
plot(lbi~as.factor(x14cgh),dat) 

#x16patcat? 
plot(lbi~as.factor(x16patcat),dat)
# New has 640 sequences, Relapse has 58, Other has 6,
#  Return .. has 2, Treatment failure has 1
# Relapse might be significantly higher than New 

#x18dot?
plot(lbi~as.factor(x18dot),dat)


plot(lbi~as.factor(x24epalcomp),dat)

plot(lbi~as.factor(x01area),dat)

plot(lbi~as.factor(l28smear),dat)

plot(lbi~as.factor(l30cultures),dat)

plot(lbi~as.factor(outcome),dat)

plot(lbi~as.factor(coughduration),dat)

plot(lbi~as.factor(tbclass),dat)

plot(lbi~as.factor(tbcategory),dat)
# return after default has 3 sequences
#  New: 636, Relapse: 62


plot(lbi~as.factor(hittbclass),dat)

plot(lbi~as.factor(peopleinhh),dat)

plot(lbi~as.factor(sleepinsameroom),dat)

plot(lbi~as.factor(smoke),dat)

plot(lbi~as.factor(relationshipwithtb),dat)

plot(lbi~as.factor(haselectricity),dat)

plot(lbi~as.factor(carmotobike),dat)

plot(lbi~as.factor(genexpertresult),dat)

plot(lbi~as.factor(rifresult),dat)

plot(lbi~as.factor(tbsymptomsidentified),dat)

plot(lbi~as.factor(lastclinicvisit),dat)

plot(lbi~as.factor(admittedhospitalone),dat)

plot(lbi~as.factor(qechwardonetbward),dat)

plot(lbi~as.factor(hospitaltimespentone),dat)

plot(lbi~as.factor(admissionunittwo),dat)

plot(lbi~as.factor(hospitaltimespenttwo),dat)

plot(lbi~as.factor(qechwardthreetbward),dat)

plot(lbi~as.factor(hospitaltimespentthree),dat)

plot(lbi~as.factor(has_worked_hs),dat)

plot(lbi~as.factor(primarycareservices),dat)

plot(lbi~as.factor(nontbward),dat)

plot(lbi~as.factor(sleepwithtb),dat)
# few datapoints in each: 23 in Yes, 21 in No


plot(lbi~as.factor(levelschool),dat)

plot(lbi~as.factor(outpatientBAH),dat)

plot(lbi~as.factor(outpatientBG),dat)

plot(lbi~as.factor(outpatientCW),dat)

plot(lbi~as.factor(outpatientMW),dat)

#hiv clinics? what do these vars mean?
plot(lbi~as.factor(hivclinicBAH),dat)

plot(lbi~as.factor(hivclinicBG),dat)

plot(lbi~as.factor(hivclinicCH),dat)

plot(lbi~as.factor(hivclinicCW),dat)

plot(lbi~as.factor(hivclinicLB),dat)
plot(lbi~as.factor(hivclinicMB),dat)

plot(lbi~as.factor(hivclinicMW),dat)

 plot(lbi~as.factor(hivclinicND),dat)

plot(lbi~as.factor(hivclinicQE),dat)

 plot(lbi~as.factor(hivclinicSL),dat)

plot(lbi~as.factor(outcome_success),dat)

#sequencing QC stuff (?)
plot(lbi~Perc_mapped_H37Rv,dat)
plot(lbi~as.factor(OriginalRun),dat)
plot(lbi~Mean_Coverage,dat)
plot(lbi~as.factor(culture_positive),dat)

