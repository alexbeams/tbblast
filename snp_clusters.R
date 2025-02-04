# find the bandwidths in LBI that maximize the signal-to-noise ratio across the posterior distribution of trees for each lineage. If you don't do this, then LBI is going to be conveying noise from uncertainty in the phylogenetic relationships, rather than a meaningful signal.

# then we find the SNP threshold for clustering that establishes the strongest correlation between LBI and
# cluster size. The clusters defined by this SNP threshold are informed by the phylogenies. 

require(Biostrings)
require(ape)
require(treeio)

source("lbi.R") #to calculate LBI

# functions that treat ambiguous sites as missing data (?)

align1 <- read.FASTA('lineage_files/alignment_lineage1.fasta')
align2 <- read.FASTA('lineage_files/alignment_lineage2.fasta')
align3 <- read.FASTA('lineage_files/alignment_lineage3.fasta')
align4 <- read.FASTA('lineage_files/alignment_lineage4.fasta')

align1.dist <- dist.dna(align1,model='N')
align2.dist <- dist.dna(align2,model='N')
align3.dist <- dist.dna(align3,model='N')
align4.dist <- dist.dna(align4,model='N')

# make clusters - use a SNP threshold of 12
align1.clusters <- cutree(hclust(align1.dist),h=12)
align2.clusters <- cutree(hclust(align2.dist),h=12)
align3.clusters <- cutree(hclust(align3.dist),h=12)
align4.clusters <- cutree(hclust(align4.dist),h=49)

getPlot <- function(cutoff,align.dist,dat){

	align.clusters <- cutree(hclust(align.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align.clusters), snpcluster = align.clusters)
	df = merge(dat, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	plot(lbi~as.double(snpclustersize),df)
}




# probably need to read in the posterior distributions of trees in here:
# refer to the bissemodels.R file

lin1phi <- read.nexus("lineage_files/beast_files/alignment1.beauti-alignment_lineage1.trees")
phi1_sample <- sample(lin1phi,size=100)

lin2phi <- read.nexus("lineage_files/beast_files/alignment2.beauti-alignment_lineage2.trees")
phi2_sample <- sample(lin2phi,size=100)

lin3phi <- read.nexus("lineage_files/beast_files/alignment3.beauti-alignment_lineage3.trees")
phi3_sample <- sample(lin3phi,size=100)



# choosing values of tau that give us clear separation among tips, and
# that tell us something slightly more useful than just distinguishing tips
# from internal nodes (as occurs when tau goes to zero):


# does it look like LBI is telling a consistent story at the tips across trees in our posterior?

#just plot lbi dist for each tip against the average lbi at each tip:

# lineage 1: two clear groups distinguished:
tau1 <- .008
align1.lbi <- sapply(phi1_sample,function(x) lbi(x,tau1))

align1.lbi.mean = apply(align1.lbi,1,mean)
plot(align1.lbi.mean[1:99],align1.lbi.mean[1:99],
	ylab='LBI from sampled trees',
	xlab='Mean LBI')
for(i in 1:100){points(align1.lbi.mean[1:99],align1.lbi[1:99,i])}
points(align1.lbi.mean[1:99],align1.lbi.mean[1:99],col='red')

mean(apply(align1.lbi,1,var))
var(apply(align1.lbi,1,mean))

var(apply(align1.lbi,1,mean)) / mean(apply(align1.lbi,1,var))

# we want to maximize the F-test statistic over tau. Each tip is a group of observations with
# as many observations as trees in the posterior sample; all groups have equal size, and
# the degrees of freedom (numerator and denominator) are constant wrt tau, so this function
# works fine for finding the optimal tau, although it isn't actually equal to the F-statistic 

getF.align1 <- function(tau){
	align1.lbi = sapply(phi1_sample,function(x) lbi(x,tau))
	align1.lbi = align1.lbi[1:99,]	
	F = var(apply(align1.lbi,1,mean)) / mean(apply(align1.lbi,1,var))
	return(F)
} 

tauvals.align1 <- seq(.005,.01,length=15)
F.align1 <- sapply(tauvals.align1, getF.align1)
plot(tauvals.align1,F.align1,type='l',xlab=bquote(tau),ylab='F',
	main='Lineage 1 LBI signal')

# lineage 2: two clear groups distinguished:
tau2 <- 5.25e-4
align2.lbi <- sapply(phi2_sample,function(x) lbi(x,tau2))

align2.lbi.mean = apply(align2.lbi,1,mean)
plot(align2.lbi.mean[1:28],align2.lbi.mean[1:28],
	xlab='Mean LBI',ylab='LBI from sampled trees')
	#ylim=c(.003,.006),xlim=c(.003,.006))
for(i in 1:100){points(align2.lbi.mean[1:28],align2.lbi[1:28,i])}
points(align2.lbi.mean[1:28],align2.lbi.mean[1:28],col='red')

mean(apply(align2.lbi,1,var))
var(apply(align2.lbi,1,mean))

var(apply(align2.lbi,1,mean)) / mean(apply(align2.lbi,1,var))

getF.align2 <- function(tau){
	align2.lbi = sapply(phi2_sample,function(x) lbi(x,tau))
	align2.lbi = align2.lbi[1:28,]	
	F = var(apply(align2.lbi,1,mean)) / mean(apply(align2.lbi,1,var))
	return(F)
} 

tauvals.align2 <- seq(.00005,.001,length=15)
F.align2 <- sapply(tauvals.align2, getF.align2)
plot(tauvals.align2,F.align2,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 2 LBI signal')


# lineage 3: two clear groups distinguished:
tau3 <- .00155
align3.lbi <- sapply(phi3_sample,function(x) lbi(x,tau3))

align3.lbi.mean = apply(align3.lbi,1,mean)
plot(align3.lbi.mean[1:61],align3.lbi.mean[1:61],
	xlab='Mean LBI',
	ylab='LBI from sampled trees')
for(i in 1:100){points(align3.lbi.mean[1:61],align3.lbi[1:61,i])}
points(align3.lbi.mean[1:61],align3.lbi.mean[1:61],col='red')

mean(apply(align3.lbi,1,var))
var(apply(align3.lbi,1,mean))

var(apply(align3.lbi,1,mean)) / mean(apply(align3.lbi,1,var))

getF.align3 <- function(tau){
	align3.lbi = sapply(phi3_sample,function(x) lbi(x,tau))
	align3.lbi = align3.lbi[1:61,]	
	F = var(apply(align3.lbi,1,mean)) / mean(apply(align3.lbi,1,var))
	return(F)
}

tauvals.align3 <- seq(.001,.003,length=15)
F.align3 <- sapply(tauvals.align3, getF.align3)
plot(tauvals.align3,F.align3,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 3 LBI signal')

# lineage 4?

lin4phi <- read.nexus('lineage_4_delphy/alex/lineage4_delphy.trees')

phi4_sample <- sample(lin4phi, 100)

# lineage 4
tau4 <- 4.7
align4.lbi <- sapply(phi4_sample,function(x) lbi(x,tau4))

align4.lbi.mean = apply(align4.lbi,1,mean)
plot(align4.lbi.mean[1:513],align4.lbi.mean[1:513],
	xlab='Mean LBI',
	ylab='LBI from sampled trees',
	main=bquote(tau[4]==.(tau4)))
for(i in 1:100){points(align4.lbi.mean[1:513],align4.lbi[1:513,i])}
points(align4.lbi.mean[1:513],align4.lbi.mean[1:513],col='red')

getF.align4 <- function(tau){
	align4.lbi = sapply(phi4_sample,function(x) lbi(x,tau))
	align4.lbi = align4.lbi[1:513,]	
	F = var(apply(align4.lbi,1,mean)) / mean(apply(align4.lbi,1,var))
	return(F)
}

tauvals.align4 <- seq(3,6,length=8)
F.align4 <- sapply(tauvals.align4, getF.align4)
plot(tauvals.align4,F.align4,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 4 LBI signal')

par(mfrow=c(2,2))
plot(tauvals.align1,F.align1,type='l',xlab=bquote(tau),ylab='F',
	main='Lineage 1 LBI signal-to-noise')

plot(tauvals.align2,F.align2,type='l',xlab=bquote(tau),ylab='F',
	main='Lineage 2 LBI signal-to-noise')

plot(tauvals.align3,F.align3,type='l',xlab=bquote(tau),ylab='F',
	main='Lineage 3 LBI signal-to-noise')

plot(tauvals.align4,F.align4,type='l',xlab=bquote(tau),ylab='F',
	main='Lineage 4 LBI signal-to-noise')

par(mfrow=c(2,2))

plot(align1.lbi.mean[1:99],align1.lbi.mean[1:99],
	ylab='LBI from sampled trees',
	xlab='Mean LBI',
	main=bquote(tau[1]==.(tau1)))
for(i in 1:100){points(align1.lbi.mean[1:99],align1.lbi[1:99,i])}
points(align1.lbi.mean[1:99],align1.lbi.mean[1:99],col='red')

plot(align2.lbi.mean[1:28],align2.lbi.mean[1:28],
	ylab='LBI from sampled trees',
	xlab='Mean LBI',
	main=bquote(tau[2]==.(tau2)))
for(i in 1:100){points(align2.lbi.mean[1:28],align2.lbi[1:28,i])}
points(align2.lbi.mean[1:28],align2.lbi.mean[1:28],col='red')


plot(align3.lbi.mean[1:61],align3.lbi.mean[1:61],
	xlab='Mean LBI',
	ylab='LBI from sampled trees',
	main=bquote(tau[3]==.(tau3)))
for(i in 1:100){points(align3.lbi.mean[1:61],align3.lbi[1:61,i])}
points(align3.lbi.mean[1:61],align3.lbi.mean[1:61],col='red')

plot(align4.lbi.mean[1:513],align4.lbi.mean[1:513],
        xlab='Mean LBI',
        ylab='LBI from sampled trees',
        main=bquote(tau[4]==.(tau4)))
for(i in 1:100){points(align4.lbi.mean[1:513],align4.lbi[1:513,i])}
points(align4.lbi.mean[1:513],align4.lbi.mean[1:513],col='red')

# need to read in dat, separate into dat1, dat2, dat3, dat4, then:

crud= data.frame(Sequence_name=names(align4.clusters), snpcluster = align4.clusters)
df = merge(dat4, crud)
df$snpclustersize = table(df$snpcluster)[df$snpcluster]

plot(lbi~as.double(snpclustersize),df)

# what SNP cutoff establishes the clearest relationship between clustersize and lbi?
# same trick as before: calculate E[ lbi | clustersize ], and choose cutoff so that
# Var[ E[lbi|clustersize]] >> E[Var[lbi|clustersize]]
# for this, we do want to be careful about correctly calculating Snedecor's F, because
# the numbers in each group of clustersize will change with the threshold
# we use the aov function this time:

getF.snpclustersize.align1 <- function(cutoff){

	align1.clusters <- cutree(hclust(align1.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align1.clusters), snpcluster = align1.clusters)
	df = merge(dat1, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	aovmod <- aov(lbi~as.double(snpclustersize),df)
	F.snpclustersize <- summary(aovmod)[[1]][["F value"]][[1]]
	p.snpclustersize <- summary(aovmod)[[1]][["Pr(>F)"]][[1]]
		
	return(c(F.snpclustersize, p.snpclustersize))

}

getF.snpclustersize.align2 <- function(cutoff){

	align2.clusters <- cutree(hclust(align2.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align2.clusters), snpcluster = align2.clusters)
	df = merge(dat2, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	aovmod <- aov(lbi~as.double(snpclustersize),df)
	F.snpclustersize <- summary(aovmod)[[1]][["F value"]][[1]]
	p.snpclustersize <- summary(aovmod)[[1]][["Pr(>F)"]][[1]]
		
	return(c(F.snpclustersize, p.snpclustersize))

}

getF.snpclustersize.align3 <- function(cutoff){

	align3.clusters <- cutree(hclust(align3.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align3.clusters), snpcluster = align3.clusters)
	df = merge(dat3, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	aovmod <- aov(lbi~as.double(snpclustersize),df)
	F.snpclustersize <- summary(aovmod)[[1]][["F value"]][[1]]
	p.snpclustersize <- summary(aovmod)[[1]][["Pr(>F)"]][[1]]
		
	return(c(F.snpclustersize, p.snpclustersize))

}


getF.snpclustersize.align4 <- function(cutoff){

	align4.clusters <- cutree(hclust(align4.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align4.clusters), snpcluster = align4.clusters)
	df = merge(dat4, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	aovmod <- aov(lbi~as.double(snpclustersize),df)
	F.snpclustersize <- summary(aovmod)[[1]][["F value"]][[1]]
	p.snpclustersize <- summary(aovmod)[[1]][["Pr(>F)"]][[1]]
		
	return(c(F.snpclustersize, p.snpclustersize))

}

cutoffs <- 10:100
F.snpclustersizes.align1 <- sapply(cutoffs, getF.snpclustersize.align1)
par(mfrow=c(2,1))
plot(cutoffs,F.snpclustersizes.align1[1,],xlab='SNP threshold', ylab='F',main='Lineage 1')
plot(cutoffs,log(F.snpclustersizes.align1[2,]),xlab='SNP threshold', ylab='log p-value')

cutoffs <- 10:100
F.snpclustersizes.align2 <- sapply(cutoffs, getF.snpclustersize.align2)
par(mfrow=c(2,1))
plot(cutoffs,F.snpclustersizes.align2[1,],xlab='SNP threshold', ylab='F',main='Lineage 2')
plot(cutoffs,log(F.snpclustersizes.align2[2,]),xlab='SNP threshold', ylab='log p-value')

cutoffs <- 10:100
F.snpclustersizes.align3 <- sapply(cutoffs, getF.snpclustersize.align3)
par(mfrow=c(2,1))
plot(cutoffs,F.snpclustersizes.align3[1,],xlab='SNP threshold', ylab='F',main='Lineage 3')
plot(cutoffs,log(F.snpclustersizes.align3[2,]),xlab='SNP threshold', ylab='log p-value')

cutoffs <- 10:100
F.snpclustersizes.align4 <- sapply(cutoffs, getF.snpclustersize.align4)
par(mfrow=c(2,1))
plot(cutoffs,F.snpclustersizes.align4[1,],xlab='SNP threshold', ylab='F',main='Lineage 4')
plot(cutoffs,log(F.snpclustersizes.align4[2,]),xlab='SNP threshold', ylab='log p-value')

# lineage 1: looks kind of pathological. LBI seems to just single out the distance sequences from the rest
# lineage 2: SNP thresholds of 15-26 establish the strongest correlation with LBI
# lineage 3: similar to lineage 1, looks pathological.  
# lineage 4: looks like any cutoff between 32 and 50, inclusive, establishes the strongest relationsip with LBI for lineage 4

#need to relabel the names of the tips in the tree according to the scheme in df (without dates)

plotLbiClust1 <- function(treeindex,cutoff){
	align1.clusters <- cutree(hclust(align1.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align1.clusters), snpcluster = align1.clusters)
	df = merge(dat1, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	nms <- sapply(df$Sequence_name, function(x) grep(x, phi1_sample[[treeindex]]$tip.label)  )
	nms <- cbind(names(nms), nms)
	nms <- nms[order(as.numeric(nms[,2])),]

	crud <- phi1_sample[[treeindex]]
	crud$tip.label <- nms[,1]

	p.tree <- ggtree(crud,layout='circular')

	p.tree <- p.tree %<+% df

	p.tree + geom_tippoint(aes(color=as.factor(snpclustersize),alpha=lbi))

}

plotLbiClust2 <- function(treeindex,cutoff){
	align2.clusters <- cutree(hclust(align2.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align2.clusters), snpcluster = align2.clusters)
	df = merge(dat2, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	nms <- sapply(df$Sequence_name, function(x) grep(x, phi2_sample[[treeindex]]$tip.label)  )
	nms <- cbind(names(nms), nms)
	nms <- nms[order(as.numeric(nms[,2])),]

	crud <- phi2_sample[[treeindex]]
	crud$tip.label <- nms[,1]

	p.tree <- ggtree(crud,layout='circular')

	p.tree <- p.tree %<+% df

	p.tree + geom_tippoint(aes(color=as.factor(snpclustersize),alpha=lbi))

}

plotLbiClust3 <- function(treeindex,cutoff){
	align3.clusters <- cutree(hclust(align3.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align3.clusters), snpcluster = align3.clusters)
	df = merge(dat3, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	nms <- sapply(df$Sequence_name, function(x) grep(x, phi3_sample[[treeindex]]$tip.label)  )
	nms <- cbind(names(nms), nms)
	nms <- nms[order(as.numeric(nms[,2])),]

	crud <- phi3_sample[[treeindex]]
	crud$tip.label <- nms[,1]

	p.tree <- ggtree(crud,layout='circular')

	p.tree <- p.tree %<+% df

	p.tree + geom_tippoint(aes(color=as.factor(snpclustersize),alpha=lbi))

}

plotLbiClust4 <- function(treeindex,cutoff){
	align4.clusters <- cutree(hclust(align4.dist),h=cutoff)
	crud= data.frame(Sequence_name=names(align4.clusters), snpcluster = align4.clusters)
	df = merge(dat4, crud)
	df$snpclustersize = table(df$snpcluster)[df$snpcluster]

	nms <- sapply(df$Sequence_name, function(x) grep(x, phi4_sample[[treeindex]]$tip.label)  )
	nms <- cbind(names(nms), nms)
	nms <- nms[order(as.numeric(nms[,2])),]

	crud <- phi4_sample[[treeindex]]
	crud$tip.label <- nms[,1]

	p.tree <- ggtree(crud,layout='circular')

	p.tree <- p.tree %<+% df

	p.tree + geom_tippoint(aes(color=as.factor(snpclustersize),alpha=lbi))

}

# what is the distribution of phylogenetic distances between sequences in the largest clusters?

plotDists <- function(snpclustersize,treeindex){
	nms <- sapply(df$Sequence_name, function(x) grep(x, phi4_sample[[treeindex]]$tip.label)  )
	nms <- cbind(names(nms), nms)
	nms <- nms[order(as.numeric(nms[,2])),]

	crud <- phi4_sample[[treeindex]]
	crud$tip.label <- nms[,1]

	# extract the pairwise distances between the sequences in the cluster
	dfcrud = df[df$snpclustersize==snpclustersize,]
	clusternms = dfcrud$Sequence_name
	
	x = cophenetic(crud)
	x=x[clusternms,clusternms]
	
	y = x[upper.tri(x)]
	hist(y)
}


