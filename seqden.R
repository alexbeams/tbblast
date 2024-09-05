rm(list=ls())

# Calculate sequence densities using patristic distances weighted by decaying exponentials
#	(similar to THD and LBI)
 
library(ape)
library(ggplot2)
library(ggtree)
library(Biostrings)
library(RColorBrewer)
library(gridExtra)
library(lme4)

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

# calculate the patristic distances
dist.tree1 <- dist.nodes(tree1)
dist.tree2 <- dist.nodes(tree2)
dist.tree3 <- dist.nodes(tree3)
dist.tree4 <- dist.nodes(tree4)

# also for the ML tree with all of the lineages:
dist.tree <- dist.nodes(tree)

dist.timetree1 <- dist.nodes(timetree1)
dist.timetree2 <- dist.nodes(timetree2)
dist.timetree3 <- dist.nodes(timetree3)
dist.timetree4 <- dist.nodes(timetree4)

# what tau values produce the most variance in the sequence density?
getseqden <- function(tree,tau){
	treedist <- dist.nodes(tree)
	tau <- tau * mean(treedist[upper.tri(treedist,diag=F)])	
	seqden <- apply(treedist,1,function(x) sum(exp(-x/tau)) )
	return( seqden)
}

#logtaus <- seq(-3,3,length=60)
#varseqdens.tree <- sapply(exp(logtaus), function(x) var(getseqden(tree,x)))
#varseqdens.tree1 <- sapply(exp(logtaus), function(x) var(getseqden(tree1,x)))
#varseqdens.tree2 <- sapply(exp(logtaus), function(x) var(getseqden(tree2,x)))
#varseqdens.tree3 <- sapply(exp(logtaus), function(x) var(getseqden(tree3,x)))
#varseqdens.tree4 <- sapply(exp(logtaus), function(x) var(getseqden(tree4,x)))
#
##rescale these for plotting
#varseqdens.tree <- varseqdens.tree/max(varseqdens.tree) 
#varseqdens.tree1 <- varseqdens.tree1/max(varseqdens.tree1) 
#varseqdens.tree2 <- varseqdens.tree2/max(varseqdens.tree2) 
#varseqdens.tree3 <- varseqdens.tree3/max(varseqdens.tree3) 
#varseqdens.tree4 <- varseqdens.tree4/max(varseqdens.tree4) 
#
#jpeg(file='untimedtree_bandwidth_variance.jpeg',
#	width=4,height=4,units='in',res=400)
#plot(logtaus,varseqdens.tree1,col='red',type='l',
#	xlab=expression(log(tau)),ylab='Scaled variance in Sequence Density',
#	main='Untimed Trees',lty=1)
#lines(logtaus,varseqdens.tree2,col='darkgreen',lty=2)
#lines(logtaus,varseqdens.tree3,col='blue',lty=3)
#lines(logtaus,varseqdens.tree4,col='purple',lty=4)
#lines(logtaus,varseqdens.tree,col='black',lty=5)
#legend('topright',legend=c(1,2,3,4,'All'),col=c('red','darkgreen','blue','purple','black'),
#	title='Lineage',lty=c(1,2,3,4,5))
#dev.off()
#
#varseqdens.timetree1 <- sapply(exp(logtaus), function(x) var(getseqden(timetree1,x)))
#varseqdens.timetree2 <- sapply(exp(logtaus), function(x) var(getseqden(timetree2,x)))
#varseqdens.timetree3 <- sapply(exp(logtaus), function(x) var(getseqden(timetree3,x)))
#varseqdens.timetree4 <- sapply(exp(logtaus), function(x) var(getseqden(timetree4,x)))
#
##rescale these for plotting
#varseqdens.timetree1 <- varseqdens.timetree1/max(varseqdens.timetree1) 
#varseqdens.timetree2 <- varseqdens.timetree2/max(varseqdens.timetree2) 
#varseqdens.timetree3 <- varseqdens.timetree3/max(varseqdens.timetree3) 
#varseqdens.timetree4 <- varseqdens.timetree4/max(varseqdens.timetree4) 
#
#jpeg(file='timedtree_bandwidth_variance.jpeg',height=4,width=4,
#	units='in',res=400)
#plot(logtaus,varseqdens.timetree1,col='red',type='l',
#	xlab=expression(log(tau)),ylab='Scaled variance in Sequence Density',
#	main='Timed Trees',lty=1)
#lines(logtaus,varseqdens.timetree2,col='darkgreen',lty=2)
#lines(logtaus,varseqdens.timetree3,col='blue',lty=3)
#lines(logtaus,varseqdens.timetree4,col='purple',lty=4)
#legend('topright',legend=c(1,2,3,4),col=c('red','darkgreen','blue','purple'),
#	title='Lineage',lty=c(1,2,3,4))
#dev.off()
#
#
#
## find the optimal bandwidths
#treeoptbw <- optimize(function(x){var(getseqden(tree,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#tree1optbw <- optimize(function(x){var(getseqden(tree1,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#tree2optbw <- optimize(function(x){var(getseqden(tree2,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#tree3optbw <- optimize(function(x){var(getseqden(tree3,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#tree4optbw <- optimize(function(x){var(getseqden(tree4,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#
#timetree1optbw <- optimize(function(x){var(getseqden(timetree1,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#timetree2optbw <- optimize(function(x){var(getseqden(timetree2,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#timetree3optbw <- optimize(function(x){var(getseqden(timetree3,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#timetree4optbw <- optimize(function(x){var(getseqden(timetree4,exp(x)))},lower=-1,upper=1,maximum=TRUE)
#
#untimed_bws <- c(tree1optbw$maximum,tree2optbw$maximum,tree3optbw$maximum,tree4optbw$maximum,treeoptbw$maximum)
#untimed_bws <- exp(untimed_bws)
#
#timed_bws <- c(timetree1optbw$maximum,timetree2optbw$maximum,timetree3optbw$maximum,timetree4optbw$maximum)
#timed_bws <- exp(timed_bws)
#
## set bandwidth parameters
#tau.tree1 <- untimed_bws[1] * mean(dist.tree1[upper.tri(dist.tree1,diag=F)])
#tau.tree2 <- untimed_bws[2] * mean(dist.tree2[upper.tri(dist.tree2,diag=F)])
#tau.tree3 <- untimed_bws[3] * mean(dist.tree3[upper.tri(dist.tree3,diag=F)])
#tau.tree4 <- untimed_bws[4] * mean(dist.tree4[upper.tri(dist.tree4,diag=F)])
#tau.tree <- untimed_bws[5] * mean(dist.tree[upper.tri(dist.tree,diag=F)])
#
#tau.timetree1 <- untimed_bws[1] * mean(dist.timetree1[upper.tri(dist.timetree1,diag=F)])
#tau.timetree2 <- untimed_bws[2] * mean(dist.timetree2[upper.tri(dist.timetree2,diag=F)])
#tau.timetree3 <- untimed_bws[3] * mean(dist.timetree3[upper.tri(dist.timetree3,diag=F)])
#tau.timetree4 <- untimed_bws[4] * mean(dist.timetree4[upper.tri(dist.timetree4,diag=F)])
#
# these bandwidths are much too large; the bandwidth giving the most variance in lineage 4 is an averge
# separation of 1000 years
tau.timetree1 <- 2
tau.timetree2 <- 2
tau.timetree3 <- 2
tau.timetree4 <- 2

tau.tree1 <- 0.0625 * mean(dist.tree1[upper.tri(dist.tree1,diag=F)])
tau.tree2 <- 0.0625 * mean(dist.tree2[upper.tri(dist.tree2,diag=F)])
tau.tree3 <- 0.0625 * mean(dist.tree3[upper.tri(dist.tree3,diag=F)])
tau.tree4 <- 0.0625 * mean(dist.tree4[upper.tri(dist.tree4,diag=F)])
tau.tree <- 0.0625 * mean(dist.tree[upper.tri(dist.tree,diag=F)])

# calculate the genome densities using exp(-d/tau)
seqden.tree <- apply(dist.tree, 1, function(x) sum(exp(-x/tau.tree)) )

seqden.tree1 <- apply(dist.tree1, 1, function(x) sum(exp(-x/tau.tree1)) )
seqden.tree2 <- apply(dist.tree2, 1, function(x) sum(exp(-x/tau.tree2)) )
seqden.tree3 <- apply(dist.tree3, 1, function(x) sum(exp(-x/tau.tree3)) )
seqden.tree4 <- apply(dist.tree4, 1, function(x) sum(exp(-x/tau.tree4)) )

seqden.timetree1 <- apply(dist.timetree1, 1, function(x) sum(exp(-x/tau.timetree1)) )
seqden.timetree2 <- apply(dist.timetree2, 1, function(x) sum(exp(-x/tau.timetree2)) )
seqden.timetree3 <- apply(dist.timetree3, 1, function(x) sum(exp(-x/tau.timetree3)) )
seqden.timetree4 <- apply(dist.timetree4, 1, function(x) sum(exp(-x/tau.timetree4)) )

#need to assign labels to the internal nodes of the trees:
tree$node.label <- paste0('node',1:tree$Nnode)
tree1$node.label <- paste0('node',1:tree1$Nnode)
tree2$node.label <- paste0('node',1:tree2$Nnode)
tree3$node.label <- paste0('node',1:tree3$Nnode)
tree4$node.label <- paste0('node',1:tree4$Nnode)

timetree1$node.label <- paste0('node',1:timetree1$Nnode)
timetree2$node.label <- paste0('node',1:timetree2$Nnode)
timetree3$node.label <- paste0('node',1:timetree3$Nnode)
timetree4$node.label <- paste0('node',1:timetree4$Nnode)

# set up dataframes with the seqden values for plotting
df.tree <- data.frame(Sequence_name = c(tree$tip.label,tree$node.label),
	seqden=seqden.tree)

df.tree1 <- data.frame(Sequence_name = c(tree1$tip.label,tree1$node.label),
	seqden=seqden.tree1)

df.tree2 <- data.frame(Sequence_name = c(tree2$tip.label,tree2$node.label),
	seqden=seqden.tree2)

df.tree3 <- data.frame(Sequence_name = c(tree3$tip.label,tree3$node.label),
	seqden=seqden.tree3)

df.tree4 <- data.frame(Sequence_name = c(tree4$tip.label,tree4$node.label),
	seqden=seqden.tree4)

df.timetree1 <- data.frame(Sequence_name = c(timetree1$tip.label,timetree1$node.label),
	seqden=seqden.timetree1)

df.timetree2 <- data.frame(Sequence_name = c(timetree2$tip.label,timetree2$node.label),
	seqden=seqden.timetree2)

df.timetree3 <- data.frame(Sequence_name = c(timetree3$tip.label,timetree3$node.label),
	seqden=seqden.timetree3)

df.timetree4 <- data.frame(Sequence_name = c(timetree4$tip.label,timetree4$node.label),
	seqden=seqden.timetree4)

# look at some other functions of distance as well, e.g. variance
df.tree$var <- apply(dist.tree, 1, var)

df.tree1$var <- apply(dist.tree1, 1, var)
df.tree2$var <- apply(dist.tree2, 1, var)
df.tree3$var <- apply(dist.tree3, 1, var)
df.tree4$var <- apply(dist.tree4, 1, var)

df.timetree1$var <- apply(dist.timetree1, 1, var)
df.timetree2$var <- apply(dist.timetree2, 1, var)
df.timetree3$var <- apply(dist.timetree3, 1, var)
df.timetree4$var <- apply(dist.timetree4, 1, var)

df.tree$mean <- apply(dist.tree, 1, mean)

df.tree1$mean <- apply(dist.tree1, 1, mean)
df.tree2$mean <- apply(dist.tree2, 1, mean)
df.tree3$mean <- apply(dist.tree3, 1, mean)
df.tree4$mean <- apply(dist.tree4, 1, mean)

df.timetree1$mean <- apply(dist.timetree1, 1, mean)
df.timetree2$mean <- apply(dist.timetree2, 1, mean)
df.timetree3$mean <- apply(dist.timetree3, 1, mean)
df.timetree4$mean <- apply(dist.timetree4, 1, mean)

# load in the metadata and the lineage/drug resistance data
drugdat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')
names(drugdat)[1] <- 'Sequence_name'

# merge them by sequence name
dat <- merge(blastdat, drugdat, by='Sequence_name')

# subset the data by lineage to generate sequence densities separetely by linage

dat1 <- dat[dat$Major.Lineage=='lineage1',]
dat2 <- dat[dat$Major.Lineage=='lineage2',]
dat3 <- dat[dat$Major.Lineage=='lineage3',]
dat4 <- dat[dat$Major.Lineage=='lineage4',]

# merge the sequence density values
datm <- merge(dat, df.tree1[,c('Sequence_name','seqden')], by='Sequence_name') 

dat1m <- merge(dat1, df.tree1[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
dat2m <- merge(dat2, df.tree2[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
dat3m <- merge(dat3, df.tree3[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
dat4m <- merge(dat4, df.timetree4[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 


# merge all of the other data onto the tree dataframes
df.tree1m <- merge(df.tree1, dat1, by='Sequence_name',all=T)
df.tree2m <- merge(df.tree2, dat2, by='Sequence_name',all=T)
df.tree3m <- merge(df.tree3, dat3, by='Sequence_name',all=T)
df.tree4m <- merge(df.tree4, dat4, by='Sequence_name',all=T)

df.timetree1m <- merge(df.timetree1, dat1, by='Sequence_name',all=T)
df.timetree2m <- merge(df.timetree2, dat2, by='Sequence_name',all=T)
df.timetree3m <- merge(df.timetree3, dat3, by='Sequence_name',all=T)
df.timetree4m <- merge(df.timetree4, dat4, by='Sequence_name',all=T)





##########################
##### plot the trees #####
##########################




# plotting parameters
treelayout <- 'circular' 
dotsize <- 1

p.tree <- ggtree(tree, layout = treelayout)

p.tree1 <- ggtree(tree1, layout = treelayout)
p.tree2 <- ggtree(tree2, layout = treelayout)
p.tree3 <- ggtree(tree3, layout = treelayout)
p.tree4 <- ggtree(tree4, layout = treelayout)

p.timetree1 <- ggtree(timetree1, layout = treelayout)
p.timetree2 <- ggtree(timetree2, layout = treelayout)
p.timetree3 <- ggtree(timetree3, layout = treelayout)
p.timetree4 <- ggtree(timetree4, layout = treelayout)

# Add the data frames to the plots
#p.tree <- p.tree %<+% df.treem

p.tree1 <- p.tree1 %<+% df.tree1m
p.tree2 <- p.tree2 %<+% df.tree2m
p.tree3 <- p.tree3 %<+% df.tree3m
p.tree4 <- p.tree4 %<+% df.tree4m

p.timetree1 <- p.timetree1 %<+% df.timetree1m
p.timetree2 <- p.timetree2 %<+% df.timetree2m
p.timetree3 <- p.timetree3 %<+% df.timetree3m
p.timetree4 <- p.timetree4 %<+% df.timetree4m


p.tree1 + geom_tippoint(aes(color=hivstatus))
p.tree4 + geom_tippoint(aes(color=hivstatus))

p.seqden.tree <- p.tree + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: ML tree')

p.seqden.tree1 <- p.tree1 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: tree 1')

p.seqden.tree2 <- p.tree2 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: tree 2')

p.seqden.tree3 <- p.tree3 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: tree 3')

p.seqden.tree4 <- p.tree4 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: tree 4')

p.seqden.timetree1 <- p.timetree1 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: timetree 1')

p.seqden.timetree2 <- p.timetree2 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: timetree 2')

p.seqden.timetree3 <- p.timetree3 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: timetree 3')

p.seqden.timetree4 <- p.timetree4 + geom_tippoint(aes(color=seqden),size=dotsize) +
	geom_nodepoint(aes(color=seqden),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Sequence density: timetree 4')

p.var.tree4 <- p.tree4 + geom_tippoint(aes(color=var),size=dotsize) +
	geom_nodepoint(aes(color=var),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Distance variance: tree 4')

p.var.timetree4 <- p.timetree4 + geom_tippoint(aes(color=var),size=dotsize) +
	geom_nodepoint(aes(color=var),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Distance variance: timetree 4')

p.mean.tree4 <- p.tree4 + geom_tippoint(aes(color=mean),size=dotsize) +
	geom_nodepoint(aes(color=mean),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Mean patristic distance: tree 4')

p.mean.timetree4 <- p.timetree4 + geom_tippoint(aes(color=mean),size=dotsize) +
	geom_nodepoint(aes(color=mean),size=dotsize) +
	scale_color_gradient(high='red',low='blue') +
	theme_tree2() +
	ggtitle('Mean patristic distance: timetree 4')

# load in the metadata and the lineage/drug resistance data
drugdat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')
names(drugdat)[1] <- 'Sequence_name'

# merge them by sequence name
dat <- merge(blastdat, drugdat, by='Sequence_name')

# subset the data by lineage to generate sequence densities separetely by linage

dat1 <- dat[dat$Major.Lineage=='lineage1',]
dat2 <- dat[dat$Major.Lineage=='lineage2',]
dat3 <- dat[dat$Major.Lineage=='lineage3',]
dat4 <- dat[dat$Major.Lineage=='lineage4',]

# append the sequence density values
dat <- merge(dat, df.tree1[,c('Sequence_name','seqden')], by='Sequence_name') 

dat1 <- merge(dat1, df.tree1[,c('Sequence_name','seqden')], by='Sequence_name') 
dat2 <- merge(dat2, df.tree2[,c('Sequence_name','seqden')], by='Sequence_name') 
dat3 <- merge(dat3, df.tree3[,c('Sequence_name','seqden')], by='Sequence_name') 
dat4 <- merge(dat4, df.timetree4[,c('Sequence_name','seqden')], by='Sequence_name') 


# create some convenient names for groups of variables

# redundant columns, unhelpful columns, 
# 	or columns with too many NA's, or minimal variation...
discardnms <- c('x01area','hittbclass','labid_duplicates','labid_january',
	'Number_paired_reads','Perc_mapped_H37Rv','Mean_Coverage',
	'OriginalRun','date','pid','tgenid','labid_links','Drug','Multi.drug')

#questions from the x01 form - filled out when clinicians register patients for the 
#	TB database
xnms <- names(dat)[5:32]

# healthcare-seeking behavior:
qechnms <- names(dat)[77:104]

# history of working in healthcare settings:
hsworknms <- names(dat)[105:113]
# only 4 people answered 'Yes'; the majority answered 'No', but lots of NAs.

# visited an outpatient clinic?
opnms <- names(dat)[116:127]

# visited an hiv clinic?
hivclinnms <- names(dat)[128:139]

# symptoms
symnms <- names(dat)[45:51]

# proxy means test for poverty:
povnms <- c('peopleinhh','sleepinsameroom','cookinglocation','smoke',
	'knowanyonewithtb','relationshipwithtb','haselectricity','fridge',
	'carmotobike','bed','radio','phone','levelschool','sleepwithtb',
	'blantyreresident')

# characteristics of the TB diagnosis
tbdxnms <- c('l28smear','l29xpert','l30cultures','l32id','smeartest','smearstatus',
	'tbclass','tbategory','genexpertresult','tbsymptomsidentified','tbcategory',
	'culture_positive')

# hiv-related characteristics
hivnms <- c('ipt','hivstatus','arvtreatment','arvduration','wereyoutakingipt')

# antibiotic information
abxnms <- c('trimoxazole','rifresult')

# demographic variables
demnms <- c('sex','age','agegroup','ageest')

# abx resistance info - second two appear redundant
resnms <- c('Drug.resistance.Tbprofiler')

# not sure what this one is - it's different from 'outcome'
huhnms <- c('outcome_success')

# the remaining columns are the sequence name, the ward ID, the outcome of treatment,
#	datecreated (differs from x05regdate) and Major.Lineage:
setdiff(names(dat), c(discardnms,xnms,qechnms,hsworknms,opnms,hivclinnms,symnms,
	povnms,tbdxnms,hivnms,abxnms,demnms,resnms,huhnms))

dat <- rbind(dat1,dat2,dat3,dat4)

# include an indicator for whether the person visited anyhivclinic, or anyop:
dat$anyhivclinic <- as.numeric(apply(dat[,hivclinnms], 1, function(x) 'Yes' %in% x))
dat$anyop <- as.numeric(apply(dat[,opnms],1,function(x) 'Yes' %in% x))


# which variables do we think matter the most? to start, aggregate all of the hivclinic
#	results into 'anyhivclinic'

modnms <- c('sex','age','hivstatus','arvtreatment','ipt','anyhivclinic','outcome',
	'tbsymptomsidentified','x04fac_code',symnms)

# some of these variables need to be reorganized before going into a regression:

# how many NAs are there for these?

sum(is.na(dat$Major.Lineage))
sum(is.na(dat$sex))
sum(is.na(dat$age))
sum(is.na(dat$hivstatus))
sum(is.na(dat$anyhivclinic))
sum(is.na(dat$outcome))
sum(is.na(dat$tbsymptomsidentified))
sum(is.na(dat$x04fac_code))

# 115 NAs in tbsymptomsidentified, 10 NAs in outcome, 
#	16 NAs in hivstatus - not too bad. 
# Do any of these variables need to be combined into a single factor?

# hivstatus and anyhivclinic are not exactly the same variable:
table(dat$hivstatus, dat$anyhivclinic)

# arvtreatment only has data points if hivstatus='Positive':
table(dat$hivstatus, dat$arvtreatment)

# coughduration obviously reported only if sympcough='Yes':
table(dat$sympcough,dat$coughduration)

# is tbsymptomsidentified a level of the responses in symnms?
dat$anysymptom <- apply(dat[,symnms],1,function(x) 'Yes' %in% x)
#	no - only 1 person in the entire dataset does not have a symptom

#how many of the symptoms listed in symnms (excluding coughduration) do people
#	in the data have?
dat$numsymptoms <- apply(dat[,symnms[1:6]],1,function(x) sum(x=='Yes'))
hist(dat$numsymptoms)
# 4 symptoms is the most common 

# 115 individuals NA for tbsymptomsidenified:
sum(is.na(dat$tbsymptomsidentified))
# sort of perplexing that this would be totally unrelated to the type of symptoms
#	reported


mod <- lm(log(seqden)~as.factor(sex) + age + as.factor(hivstatus) +
	as.factor(anyhivclinic) + as.factor(outcome) + as.factor(tbsymptomsidentified) + 
	as.factor(x04fac_code) + numsymptoms, dat4 )

# look at the variables in isolation for lineage 4:

# sex and age might matter on their own:
summary(lm(log(seqden)~as.factor(sex),dat4))

# this may be meaningless. The magnitude of the effect is pretty small
summary(lm(log(seqden)~age,dat4))

# no signal from hivstatus:
summary(lm(log(seqden)~as.factor(hivstatus),dat4))

# no signal from anyhivclinic:
summary(lm(log(seqden)~as.factor(anyhivclinic),dat4))

# no signal from outcome:
summary(lm(log(seqden)~as.factor(outcome),dat4))

# no signal from tbsymptomsidentified (but this had 115 NAs):
summary(lm(log(seqden)~as.factor(tbsymptomsidentified),dat4))

# GateWay from x04fac_code is where the most immediately related sequences comes from (using 5yr bandwidth)
summary(lm(log(seqden)~as.factor(x04fac_code),dat4))

# strong signal from Drug.resistance.Tbprofiler, indicates that the highest density sequences
#	are Sensitive (only at larger bw's, however):
summary(lm(log(seqden)~as.factor(Drug.resistance.Tbprofiler),dat4))
# However: of 701 rows, 667 sensitivee, 30 resistant, 4 are MDR
table(dat4$Drug.resistance.Tbprofiler)

# no signal from anyop:
summary(lm(log(seqden)~as.factor(anyop),dat4))


# relapse in tbcategory seems positively associated with seqden at 5-yr bw:
summary(lm(log(seqden)~as.factor(tbcategory),dat4))

# no signal from tbclass
summary(lm(log(seqden)~as.factor(tbclass),dat4))

# no signal from l28smear 
summary(lm(log(seqden)~as.factor(l28smear),dat4))


# what about lineage 1?

# sex and age might matter on their own:
summary(lm(log(seqden)~as.factor(sex),dat1))

# this may be meaningless. The magnitude of the effect is pretty small
summary(lm(log(seqden)~age,dat1))

# no signal from hivstatus:
summary(lm(log(seqden)~as.factor(hivstatus),dat1))

# lineage 4 has a significantly smaller share of HIV:
summary(glm(as.factor(hivstatus)~as.factor(Major.Lineage),family=binomial(link='logit'),data=dat))

# there does appear to a signal with hivstatus over all lineages:
summary(lm(log(seqden)~as.factor(hivstatus),dat))

# removing just lineage 4 doesn't reveal a relationshp with hivstatus:
summary(lm(log(seqden)~as.factor(hivstatus),subset(dat, Major.Lineage %in% c('lineage1','lineage2','lineage3') ) ))



# no signal from anyhivclinic:
summary(lm(log(seqden)~as.factor(anyhivclinic),dat1))

# no signal from outcome:
summary(lm(log(seqden)~as.factor(outcome),dat1))

# no signal from tbsymptomsidentified (but this had 115 NAs):
summary(lm(log(seqden)~as.factor(tbsymptomsidentified),dat1))

#  signals x04fac_code at 5-yr bw: BT Adventist has low density, Chilomoni, Queen Elizabeth, Zingwangwa have higher density 
summary(lm(log(seqden)~as.factor(x04fac_code),dat1))

# strong signal from Drug.resistance.Tbprofiler, indicates that the highest density sequences
#	are Sensitive:
summary(lm(log(seqden)~as.factor(Drug.resistance.Tbprofiler),dat1))
# However: of 99 rows, 91 sensitivee, 8 resistant, 0 are MDR
table(dat1$Drug.resistance.Tbprofiler)

# no signal from anyop:
summary(lm(log(seqden)~as.factor(anyop),dat1))


# no signal from tbcategory:
summary(lm(log(seqden)~as.factor(tbcategory),dat1))

# no signal from tbclass
summary(lm(log(seqden)~as.factor(tbclass),dat1))

# no signal from l28smear 
summary(lm(log(seqden)~as.factor(l28smear),dat1))


# what about numsymptoms?
summary(lm(log(seqden)~numsymptoms,dat2))


# some other interesting plots:

plot(seqden~as.factor(smearstatus),dat4)

plot(seqden~as.factor(Major.Lineage),dat)

plot(seqden~as.factor(Drug.resistance.Tbprofiler),dat)

plot(seqden~as.factor(outcome),dat)

plot(seqden~as.factor(agegroup),dat)

plot(seqden~as.factor(hivstatus),dat)

plot(seqden~as.factor(hivstatus),dat2)


plot(seqden~as.factor(ipt),dat3)
plot(seqden~as.factor(ipt),dat)

# Unknown_Inconclusive only has 2 entries
plot(seqden~as.factor(smearstatus),dat)

# are the relationships different for different lineages?
plot(seqden~as.factor(trimoxazole),dat)

plot(seqden~as.factor(rifresult),dat1)
plot(seqden~as.factor(rifresult),dat4)



