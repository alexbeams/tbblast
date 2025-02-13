rm(list=ls())

# use LBI and clustering to examine relationships between recent transmission
# success and drug resistance.
# Show that using statistics like LBI is a more powerful approach than clustering.

require(ape)
require(ggtree)

# read in the pre-packaged lbi and thd functions for comparison with the 
# simpler exp(-x/tau)

source('lbi.R')
source('thd.R')

# read in the timed trees for the Blantyre data produced from IQ-tree:

blantyre.timetree1 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage1_results/timetree.nexus')
blantyre.timetree2 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage2_results/timetree.nexus')
blantyre.timetree3 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage3_results/timetree.nexus')
blantyre.timetree4 <- ape::read.nexus('constant_site_correction_kimura/treetime_lineage4_results/timetree.nexus')

# also, let's read in the posterior samples for those same lineages obtained
# from BEAST (lineages 1-3) or Delphy (lineage 4):

blantyre.lin1phi <- read.nexus("lineage_files/beast_files/alignment1.beauti-alignment_lineage1.trees")
blantyre.phi1_sample <- sample(blantyre.lin1phi,size=100)

blantyre.lin2phi <- read.nexus("lineage_files/beast_files/alignment2.beauti-alignment_lineage2.trees")
blantyre.phi2sample <- sample(blantyre.lin2phi,size=100)

blantyre.lin3phi <- read.nexus("lineage_files/beast_files/alignment3.beauti-alignment_lineage3.trees")
blantyre.phi3_sample <- sample(blantyre.lin3phi,size=100)

blantyre.lin4phi <- read.nexus('lineage_4_delphy/alex/lineage4_delphy.trees')
blantyre.phi4_sample <- sample(blantyre.lin4phi, 100)

# read in the alignments for each lineage to make SNP clusters:
blantyre.align1 <- read.FASTA('lineage_files/alignment_lineage1.fasta')
blantyre.align2 <- read.FASTA('lineage_files/alignment_lineage2.fasta')
blantyre.align3 <- read.FASTA('lineage_files/alignment_lineage3.fasta')
blantyre.align4 <- read.FASTA('lineage_files/alignment_lineage4.fasta')

blantyre.align1.dist <- dist.dna(blantyre.align1,model='N')
blantyre.align2.dist <- dist.dna(blantyre.align2,model='N')
blantyre.align3.dist <- dist.dna(blantyre.align3,model='N')
blantyre.align4.dist <- dist.dna(blantyre.align4,model='N')

# read in the metadata for the Blantyre sequences:

# load in the metadata and the lineage/drug resistance data
blantyre.drugdat <- read.csv('Malawi_final_stats.csv')
blantyre.blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')
names(blantyre.drugdat)[1] <- 'Sequence_name'

# merge them by sequence name
blantyre.dat <- merge(blantyre.blastdat, blantyre.drugdat, by='Sequence_name')

# load in the dataframe that contains "cluster", a nbd identifier
require(readxl)
blast <- read_excel('BLAST.xls')

# each of dat and blast has pid's that are not in the other (551 in common)
blantyre.dat <- merge(blantyre.dat, blast[,c('pid','cluster')], all.x=T)

# subset the data by lineage to generate sequence densities separetely by linage

blantyre.dat1 <- blantyre.dat[blantyre.dat$Major.Lineage=='lineage1',]
blantyre.dat2 <- blantyre.dat[blantyre.dat$Major.Lineage=='lineage2',]
blantyre.dat3 <- blantyre.dat[blantyre.dat$Major.Lineage=='lineage3',]
blantyre.dat4 <- blantyre.dat[blantyre.dat$Major.Lineage=='lineage4',]

# now, read in the tree from Moldova and its metadata:

moldova.tree <- read.tree("/Users/abeams/Documents/Moldova_TB/Moldova_pure_noAMR.fasta.treefile")
moldova.dat <- read.csv("/Users/abeams/Documents/Moldova_TB/Moldova_pure_metadata.csv")


# calculate bandwidths, tau, in two different ways:
# 1. choose the tau value that establishes the strongest positive correlation between
#	SNP cluster size and LBI for a given snpthreshold
# 2. choose the tau that maximizes the signal-to-noise ratio for LBI across posterior samples of trees.
#	Q: what SNP threshold will this correspond to?
# 3. choose the tau that maximizes the signal-to-noise ratio for LBI across levels of a given treatment
# 4. Problem: want to identify the SNP threshold that the data have the most information about. 
#	However, fixing tau and changing snpthreshold gives F-statistics that cannot be compared with
#	each other. Think the best we can do is choose some grid of snpthresholds, find the tau's
#	that appear to correspond to those, and compare those tau's to the ones that the posterior
#	sample of trees supports. Then we can say that the data are powered to detect differences
#	between clusters with thresholds of x

# let's start with the Moldova data first (no posterior available currently)

# calculate SNP clusters
moldova.align <- read.FASTA('/Users/abeams/Documents/Moldova_TB/Moldova_Pure2067_AMR_SNPs.fasta')
moldova.align.dist <- dist.dna(moldova.align,model='N')
moldova.align.clusters <- cutree(hclust(moldova.align.dist),h=12) # <- SNP threshold here

# the lbi function returns NA's on moldova.tree, probably because it's unrooted. Let's worry 
# about that later, and just use exp(-x/tau) for now

getSeqden <- function(tree, tau){
	x <- cophenetic(tree)
	seqden <- apply(x,1,function(x) sum(exp(-x/tau)))
	return(seqden)
}

getClusts <- function(align.dist,snpthresh){
	clusters = cutree(hclust(align.dist),h=snpthresh)
	return(clusters)
}

getSeqdenvsClusts <- function(tree,tau,align.dist,snpthresh){
	seqden = getSeqden(tree,tau)
	clusters = getClusts(align.dist,snpthresh)
	seqden = data.frame(SRR = names(seqden), seqden = seqden)
	clusters = data.frame(SRR = names(clusters), clusters=clusters)
	dat = merge(seqden,clusters)
	dat$snpclustersize = table(dat$clusters)[dat$clusters]	
	return(dat)
}

# if we use 0.0625 mean(cophenetic(moldova.tree)), what SNP threshold is this simlar to?
# ! This move is illegal
# use aov to calculate an F-statistic as a functino of snpthreshold
getF <- function(snpthresh){
	x = getSeqdenvsClusts(moldova.tree,0.0625*mean(cophenetic(moldova.tree)),moldova.align.dist,snpthresh)
	f = summary(aov(seqden~snpclustersize,x))[[1]][1,'F value']
	return(f)
}

snpthreshvals <- c(40:100,by=5)
fvals <- sapply(snpthreshvals, getF)

par(mfrow=c(2,1))
plot(snpthreshvals,fvals,main='F-statistic',xlab='SNP threshold',ylab='F-statistic')
x = getSeqdenvsClusts(moldova.tree,0.0625 *mean(cophenetic(moldova.tree)),moldova.align.dist,60)
boxplot(seqden~snpclustersize,x,main='SNP threshold = 60')


# Work the other way: what tau corresponds to a SNP threshold of 12? of 5?

getF <- function(tau){
	x = getSeqdenvsClusts(moldova.tree,tau,moldova.align.dist,12)
	f = summary(aov(seqden~snpclustersize,x))[[1]][1,'F value']
	return(f)
}

logtauvals <- seq(-12,-9,length=60)
fvals <- sapply(exp(logtauvals),getF)
par(mfrow=c(2,1))
plot(logtauvals,fvals,xlab=bquote(log(tau)), ylab='F-statistic')
x = getSeqdenvsClusts(moldova.tree, exp(-10.75), moldova.align.dist,12)
boxplot(seqden~snpclustersize,x)

# can also work the other way: what tau corresponds to a SNP threshold of 12? of 5?

getF <- function(tau){
	x = getSeqdenvsClusts(moldova.tree,tau,moldova.align.dist,5)
	f = summary(aov(seqden~snpclustersize,x))[[1]][1,'F value']
	return(f)
}

logtauvals <- seq(-12,-9,length=20)
fvals <- sapply(exp(logtauvals),getF)
par(mfrow=c(2,2))
plot(logtauvals,fvals,xlab=bquote(log(tau)), ylab='F-statistic')
x = getSeqdenvsClusts(moldova.tree, exp(-11.35), moldova.align.dist,5)
boxplot(seqden~snpclustersize,x,main='log tau = -11.35 / SNP threshold = 5')

crud = merge(moldova.dat,x)

boxplot(snpclustersize~DR_status,crud)
boxplot(seqden~DR_status,crud)

###
# Result 1: #
###

# at a 5-SNP threshold, snpclustersize~DR_status is barely significant:
kruskal.test(snpclustersize~DR_status,crud)
# pairwise comparisions indicate Sensitive is significantly different from MDR_TB, but no others:
pairwise.wilcox.test(crud$snpclustersize,crud$DR_status,method.adjust='bonferroni')

# At the corresponding tau-value, however, seqden is unquestionably different across
# the different DR_status levels, with MDR-TB clearly having higher seqden than Sensitive
kruskal.test(seqden~DR_status,crud)
pairwise.wilcox.test(crud$seqden,crud$DR_status,method.adjust='bonferroni')
# Sensitive is significantly different from all groups except RR-TB and XDR-TB (small numbers of those)

# is there a particular tau-value that establishes the most obvious difference across DR_status?

getF <- function(tau){
	x = getSeqdenvsClusts(moldova.tree,tau,moldova.align.dist,5)
	dat = merge(moldova.dat,x)	
	f = summary(aov(seqden~DR_status,dat))[[1]][1,'F value']
	return(f)
}

logtauvals <- seq(-10,-7,length=20)
fvals <- sapply(exp(logtauvals), getF)
plot(logtauvals,fvals)

par(mfrow=c(2,1))
plot(logtauvals,fvals,xlab=bquote(log(tau)), ylab='F-statistic')
x = getSeqdenvsClusts(moldova.tree, exp(-8.5), moldova.align.dist,5)
crud = merge(moldova.dat,x)
boxplot(seqden~DR_status,crud,main='log tau = -8.5 ')

# what SNP threshold does this correspond to?
# ! This is probably not valid to do; as we change snpthreshold, the number of points in each
# snpclustersize will change, so that the different F-statistics aren't really comparable

getF <- function(snpthreshold){
	x = getSeqdenvsClusts(moldova.tree,exp(-8.5),moldova.align.dist,snpthreshold)
	dat = merge(moldova.dat,x)	
	f = summary(aov(seqden~snpclustersize,dat))[[1]][1,'F value']
	return(f)
}


snpthresholdvals <- seq(40,150,by=5)
fvals <- sapply(snpthresholdvals, getF)
par(mfrow=c(2,1))
plot(snpthresholdvals,fvals,main='exp(-8.5)',xlab='SNP cluster threshold',ylab='F-statistic: seqden~snpclusterize')
x = getSeqdenvsClusts(moldova.tree, exp(-8.5), moldova.align.dist,60)
crud=merge(moldova.dat,x)
boxplot(seqden~snpclustersize,crud,main='SNP threshold=60')

# tau=exp(-8.5) maximizes the association between seqden and DR_status; this corresponds to a SNP-threshold
# of 50-80 or so. With this threshold, there is a pretty clear correlation between seqden and snpclustersize,
# but it isn't perhaps as clean as it could be

# Drug resistance is actively growing, because MDR-TB is more commonly found
# in larger SNP clusters at thresholds of 12. The virtue of using an LBI-type statistic is that this
# still holds at smaller bandwidths commensurate with clusters of 5-SNP thresholds, where an association of
# snpclustersize and DR_status starts to break down. The data have the most power distinguish seqden accordin 
# to DR_status at larger bandwidths corresponding to clusters with SNP-thresholds of 50-80, which is
# very similar to the value used in the Moldova paper (0.0625 * mean(cophenetic(moldova.tree))).
# An LBI-type statistic always appears to be more powerful to distinguish differences among DR_status
# compared to snpclustersize, and this is particularly useful at low bandwidths that correspond to 
# recent transmission.

# Just for fun: what happens if we sub-sample the Moldova lineage 4 datato be the same size as the Blantyre lineage 4
# data, and we also keep the same number of resistant sequences?

moldova.dat4 = moldova.dat[moldova.dat$Lineage=='lineage4',]

#remove the "Other" for DR_status:
moldova.dat4 = moldova.dat4[moldova.dat4$DR_status != 'Other',]

getP <- function(){
	moldova.dat4.sensitive = moldova.dat4[moldova.dat4$DR_status=='Sensitive',]
	moldova.dat4.resistant = moldova.dat4[moldova.dat4$DR_status!='Sensitive',]

	# now, sub-sample these to match the numbers in Blantyre:
	crud.moldova.sensitive = moldova.dat4.sensitive[sample(687,496),]
	crud.moldova.resistant = moldova.dat4.resistant[sample(419,17),]
	# just name all of the resitant categories "Resistant"
	crud.moldova.resistant$DR_status <- "Resistant"

	crud.moldova.dat4 = rbind(crud.moldova.sensitive,crud.moldova.resistant)
	sampnms = crud.moldova.dat4[,'SRR']

	y = x[sampnms,sampnms]

	seqden.moldova.samp = apply(y,1,function(x) sum(exp(-x/.0001)))
	seqden.moldova.samp=data.frame(SRR=names(seqden.moldova.samp),seqden=seqden.moldova.samp)
	cruddat = merge(crud.moldova.dat4,seqden.moldova.samp)
	boxplot(seqden~DR_status,cruddat)
	p = wilcox.test(seqden~DR_status,cruddat)[[3]]
	return(p)
}

# simulating this 1000 times shows that about 17% of the p-values fall below .05; so there is still a tiny
# bit of a signal from resistance, even with this extreme sub-sampling We have of course lost
# most of the power to detect a difference. However, the fact that we reject the null hypothesis
# more than 5% of the time indicates that phylogenetic relationships independent of the number of observed
# sequences with resistance do convey some information about the benefit of drug resistance
# although, must admit that we are employing a 2-sided test so Resistant might not always have higher seqden
# when we reject the distributions being equal
# Practially, however, the lack of resistant sequences means that the vast majority of the time we will not
# be rejcting the null hypothesis even if we should


# What do the Blantyre data show?

#first, we need to extract the lin4nms from the data, because blantyre.phi4_sample has the Delphy
# naming convention with dates

blantyre.lin4nms = blantyre.dat4$Sequence_name
y <- sapply(blantyre.phi4_sample[[1]]$tip.label, function(x)  grep(x, blantyre.lin4nms)  ) 

tree <- blantyre.phi4_sample[[1]]
tree$tip.label = blantyre.lin4nms[y]

x = cophenetic(tree)

