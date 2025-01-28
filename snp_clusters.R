# find the bandwidths in LBI that maximize the signal-to-noise ratio across the posterior distribution of trees for each lineage.


require(Biostrings)
require(ape)
require(treeio)

source("lbi.R") #to calculate LBI

##load in the alignments - these are just the variable sites:
#align1 <- readDNAStringSet('lineage_files/alignment_lineage1.fasta')
#align2 <- readDNAStringSet('lineage_files/alignment_lineage2.fasta')
#align3 <- readDNAStringSet('lineage_files/alignment_lineage3.fasta')
#align4 <- readDNAStringSet('lineage_files/alignment_lineage4.fasta')
#
## compare with the full alignment - check that we didn't mess up the lineage
## files when we split them up into different lineages:
#align.full <- readDNAStringSet('Malawi_final_filtered.fasta')
#
##convert these to matrices:
#align1mat <- as.matrix(align1)
#align2mat <- as.matrix(align2)
#align3mat <- as.matrix(align3)
#align4mat <- as.matrix(align4)
#
## calculate the pairwise SNP distances for alignment 1:
#align1dists <- matrix(nrow=nrow(align1mat),ncol=nrow(align1mat))
#
## just write a dumb nested for loop to fill in the matrix (it's symmetric): 
#for(i in 1:nrow(align1mat)){
#	for(j in i:nrow(align1mat)){
#		align1dists[i,j] <- sum(align1mat[i,] != align1mat[j,])
#		align1dists[j,i] <- align1dists[i,j]	
#	}
#}
#
## alignment 2:
#align2dists <- matrix(nrow=nrow(align2mat),ncol=nrow(align2mat))
#
## just write a dumb nested for loop to fill in the matrix (it's symmetric): 
#for(i in 1:nrow(align2mat)){
#	for(j in i:nrow(align2mat)){
#		align2dists[i,j] <- sum(align2mat[i,] != align2mat[j,])
#		align2dists[j,i] <- align2dists[i,j]	
#	}
#}
#
#
#align3dists <- matrix(nrow=nrow(align3mat),ncol=nrow(align3mat))
#
## just write a dumb nested for loop to fill in the matrix (it's symmetric): 
#for(i in 1:nrow(align3mat)){
#	for(j in i:nrow(align3mat)){
#		align3dists[i,j] <- sum(align3mat[i,] != align3mat[j,])
#		align3dists[j,i] <- align3dists[i,j]	
#	}
#}
#
#align4dists <- matrix(nrow=nrow(align4mat),ncol=nrow(align4mat))
#
## just write a dumb nested for loop to fill in the matrix (it's symmetric): 
#for(i in 1:nrow(align4mat)){
#	for(j in i:nrow(align4mat)){
#		align4dists[i,j] <- sum(align4mat[i,] != align4mat[j,])
#		align4dists[j,i] <- align4dists[i,j]	
#	}
#}
#
## compare with the function built in to Biostrings (appears to match, but this one runs much faster):
#align1dists.compare <-  pwalign::stringDist(align1,method='hamming')
#align2dists.compare <-  pwalign::stringDist(align2,method='hamming')
#align3dists.compare <-  pwalign::stringDist(align3,method='hamming')
#align4dists.compare <-  pwalign::stringDist(align4,method='hamming')
#
## look at the full alignment
#aligndists.compare <- pwalign::stringDist(align.full,method='hamming')
# 
## the minimum hamming distances look pretty large:
#min(align1dists.compare)
#min(align2dists.compare)
#min(align3dists.compare)
#min(align4dists.compare)
#
## indeed, the observed min of 31 is also true for the original alignment
#min(aligndists.compare)
#
## what about if we examine sequences collected from the same individuals at different times?
#align.allsequences <- readDNAStringSet('Malawi_final.fasta')
#
#align.withinhost <- setdiff(align.allsequences,align.full)
#align.withinhost.compare <- pwalign::stringDist(align.withinhost,method='hamming')
#
## distribution of distances looks interesting, esp. for lineage 4:
#par(mfrow=c(2,2))
#hist(align1dists.compare,breaks=50)
#hist(align2dists.compare)
#hist(align3dists.compare)
#hist(align4dists.compare,breaks=100)
#
## so it looks like there are 3 distinct groups within lineage 4, 2 distinct groups within lineage 1 
#
#heatmap(as.matrix(align1dists.compare), symm=T)
#heatmap(as.matrix(align2dists.compare), symm=T)
#heatmap(as.matrix(align3dists.compare), symm=T)
#heatmap(as.matrix(align4dists.compare), symm=T)

# everything prior to this treats ambiguous sites as nucleotides; let's use the ape 
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
align4.clusters <- cutree(hclust(align4.dist),h=12)

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

getF.align1 <- function(tau){
	align1.lbi = sapply(phi1_sample,function(x) lbi(x,tau))
	align1.lbi = align1.lbi[1:99,]	
	F = var(apply(align1.lbi,1,mean)) / mean(apply(align1.lbi,1,var))
	return(F)
} 

tauvals.align1 <- seq(.005,.01,length=15)
F.align1 <- sapply(tauvals.align1, getF.align1)
plot(tauvals.align1,F.align1,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 1 LBI signal')

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
	main='Alignment 1 LBI signal-to-noise')

plot(tauvals.align2,F.align2,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 2 LBI signal-to-noise')

plot(tauvals.align3,F.align3,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 3 LBI signal-to-noise')

plot(tauvals.align4,F.align4,type='l',xlab=bquote(tau),ylab='F',
	main='Alignment 4 LBI signal-to-noise')

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


