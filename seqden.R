rm(list=ls())

# Calculate sequence densities using patristic distances weighted by decaying exponentials
#	(similar to THD and LBI)
 
library(ape)
library(ggplot2)
library(ggtree)
#library(Biostrings)
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
tau.timetree1 <- 5
tau.timetree2 <- 5
tau.timetree3 <- 5
tau.timetree4 <- 5

# no idea what the bandwidths for untimed trees should be; probably don't even want to use
# the untimed trees
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
	'culture_negitive')

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


# include an indicator for whether the person visited anyhivclinic, or anyop:
dat$anyhivclinic <- as.numeric(apply(dat[,hivclinnms], 1, function(x) 'Yes' %in% x))
dat$anyop <- as.numeric(apply(dat[,opnms],1,function(x) 'Yes' %in% x))

#how many NAs?
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


# subset the data by lineage to generate sequence densities separetely by linage

dat1 <- dat[dat$Major.Lineage=='lineage1',]
dat2 <- dat[dat$Major.Lineage=='lineage2',]
dat3 <- dat[dat$Major.Lineage=='lineage3',]
dat4 <- dat[dat$Major.Lineage=='lineage4',]

# merge the sequence density values - these now have a row for each internal node, all with NA except
# for the seqden values

dat1m <- merge(dat1, df.timetree1[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
dat2m <- merge(dat2, df.timetree2[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
dat3m <- merge(dat3, df.timetree3[,c('Sequence_name','seqden')], by='Sequence_name',all=T) 
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

# now reproduce dat with all lineages in it to do regressions later, and might as well
# leave out all of the rows for internal nodes because those don't have any information
dat1 <- dat1m[dat1m$Sequence_name %in% dat$Sequence_name,] 
dat2 <- dat2m[dat2m$Sequence_name %in% dat$Sequence_name,] 
dat3 <- dat3m[dat3m$Sequence_name %in% dat$Sequence_name,] 
dat4 <- dat4m[dat4m$Sequence_name %in% dat$Sequence_name,] 

datold = dat
dat = rbind(dat1,dat2,dat3,dat4)


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

###############
# Regressions #
###############


modnms <- c('sex','age','hivstatus','arvtreatment','ipt','anyhivclinic','outcome',
	'tbsymptomsidentified','x04fac_code',symnms)


mod <- lm(log(seqden)~as.factor(sex) + age + as.factor(hivstatus) +
	as.factor(anyhivclinic) + as.factor(outcome) + as.factor(tbsymptomsidentified) + 
	as.factor(x04fac_code) + numsymptoms, dat4 )

# look at the variables in isolation for lineage 4:

# sex and age might matter on their own at larger bandwidths, but not for small ones:
summary(lm(log(seqden)~as.factor(sex),dat4))

# this may be meaningless. The magnitude of the effect is pretty small (again, only for a larger bandwidth)
summary(lm(log(seqden)~age,dat4))

# no signal from hivstatus:
summary(lm(log(seqden)~as.factor(hivstatus),dat4))

# no signal from anyhivclinic:
summary(lm(log(seqden)~as.factor(anyhivclinic),dat4))

# no signal from outcome:
summary(lm(log(seqden)~as.factor(outcome),dat4))

# no signal from tbsymptomsidentified (but this had 115 NAs):
summary(lm(log(seqden)~as.factor(tbsymptomsidentified),dat4))

# GateWay from x04fac_code is where the most immediately related sequences comes from (using 2yr bandwidth)
summary(lm(log(seqden)~as.factor(x04fac_code),dat4))

# strong signal from Drug.resistance.Tbprofiler, indicates that the highest density sequences
#	are Sensitive (only at larger bw's, however):
summary(lm(log(seqden)~as.factor(Drug.resistance.Tbprofiler),dat4))
# However: of 701 rows, 667 sensitivee, 30 resistant, 4 are MDR
table(dat4$Drug.resistance.Tbprofiler)

# no signal from anyop:
summary(lm(log(seqden)~as.factor(anyop),dat4))


# relapse in tbcategory seems positively associated with seqden at 2-yr bw:
summary(lm(log(seqden)~as.factor(tbcategory),dat4))
summary(lm(log(seqden)~as.factor(tbcategory)+Major.Lineage,dat))
# relationship is clearer including all lineages; the effect decreases a small amount

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
# at 2-yr bw, Chilomoni and GatWay have significantly higher density for lineage 2; for lineage 4, Chilomoni, Gateway,
#	Ndirande, and Queen Elizabeth Central Hospital have significantly elevated density; no signal for lineages 1 or 3

summary(lm(log(seqden)~as.factor(x04fac_code)+Major.Lineage,dat))
# at 2-yr bw, Chilomoni, Gateway, Queen Elizabeth Central are almost signifigant (higher density)

# strong signal from Drug.resistance.Tbprofiler, indicates that the highest density sequences
#	are Sensitive (but only for lineage 1):
summary(lm(log(seqden)~as.factor(Drug.resistance.Tbprofiler),dat1))
# However: of 99 rows, 91 sensitivee, 8 resistant, 0 are MDR
table(dat1$Drug.resistance.Tbprofiler)

summary(lm(log(seqden)~as.factor(Drug.resistance.Tbprofiler)+Major.Lineage,dat))
# no clear signal here

# no signal from anyop:
summary(lm(log(seqden)~as.factor(anyop),dat1))


# weak signal from tbcategory (2-yr bw):
summary(lm(log(seqden)~as.factor(tbcategory),dat1))

summary(lm(log(seqden)~as.factor(tbcategory)+Major.Lineage,dat))
# strongest signal from looking at the full data (lineages 1 and 4 show a positive relationship with seqden and relapse
#	on their own)

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


summary(lm(log(seqden)~as.factor(tbcategory)+Major.Lineage+x04fac_code+age+sex+l28smear,dat))

relapsemod <- lm(log(seqden)~tbcategory+Major.Lineage+x04fac_code+age+sex+l28smear+hivstatus+anyhivclinic+outcome+tbsymptomsidentified+sympsweat+sympfever+sympblood+sympcough+ipt+wereyoutakingipt,dat)


###################################
#### GLMs for sequence desnity ####
###################################


# let's start with variables that don't have any missingness
missingness <- apply(dat,2,function(x) sum(is.na(x)))
completes <- names(dat)[which(missingness==0)]

# some redundancies in completes, but this is the big model to start from:
mod <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))

par(mfrow=c(2,2))
# residual plots don't look fantastic, but certainly better than using a linear model for log(seqden). QQ plot isn't expected to fall
# on the 1-1 line since the error structure is no longer Gaussian
plot(mod)

# different corrections like bonferroni, fdr, etc. identify lineage as important, but nothing else:
which(p.adjust(coef(summary(mod))[,4], method='hochberg')<.05)
which(p.adjust(coef(summary(mod))[,4], method='bonferroni')<.05)
which(p.adjust(coef(summary(mod))[,4], method='fdr')<.05)
which(p.adjust(coef(summary(mod))[,4], method='BY')<.05)


# try removing each variable (or groups of variables, such as symptoms and poverty) one at a time to see how much the fit deteriorates
# to identify the most important contributions to the model fit
# use AIC and/or ANOVA to assess model performance

# x04fac_code removed:
mod1 <- glm(seqden~x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod, mod1)
anova(mod,mod1)
# no major change

# x05year removed (wasn't significant on its own in the summary, either):
mod1 <- glm(seqden~x04fac_code           +x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod, mod1)
# no major change

# x09iptpst removed:
mod1 <- glm(seqden~x04fac_code+x05year+             x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# no major change

# x10iptdiag removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+              x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x10iptdiag can be removed

# x17tbtype removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+              x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x17tbtype can be removed

# x16patcat removed
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+              x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x16patcat can be removed

# x20hiv removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+           x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x20hiv can be removed

# x21art removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+     x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x21art can be removed

# x22cotri removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+        x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x22cotri can be removed

# x02res removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+      sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# x02res can be removed

# sex removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+   age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# sex can be removed

# age removed (and wasn't significant on its own in the summary):
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+    smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# age can be removed

# smeartest removed:
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+          sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# smeartest can be removed

# try removing the symptom variables
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+   fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
# removing symptoms variables makes AIC increase by 9 points; anova suggests symptoms aren't that important
AIC(mod,mod1)
anova(mod,mod1)


# try removing economic variables to see how the fit compares
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+   Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
#anova agrees:
anova(mod,mod1)
# AIC and anova suggest the economic variables are important on their own

# try removing Major.Lineage
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+        culture_positive+Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# can't remove Major.Lineage

# try removing culture_positive
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+        Drug.resistance.Tbprofiler,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# can remove culture_positive (but tried changing bandwidths to 2 years; now this seems important) 

# try removing Drug.resistance.Tbprofiler
mod1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# can remove culture_positive 

# Lots of these made borderline important contributions to model fit. Major.Lineage and the economic variables were the most important, 
# i.e. only ones which make an obvious contribution by themselves
# First, ran this for seqden with a 5-yr bandwidth. Running at 2-yr bandwidth, culture_positive seems to become more important, and
# the result for Major.Lineage and economic variables still holds

# what if we compare the full model with just the ones with economic variables and lineage?
modlineco <- glm(seqden~fridge+carmotobike+bed+radio+phone+Major.Lineage,dat,family=Gamma(link='identity'))
AIC(mod,modlineco)
anova(mod,modlineco)
# it seems to be doing pretty well

# what about just Major.Lineage?
mod1 <- glm(seqden~Major.Lineage,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# very pronounced difference - indicates we need the economic data to help the model fit

# what about just the economic data?
mod1 <- glm(seqden~fridge+carmotobike+bed+radio+phone,dat,family=Gamma(link='identity'))
AIC(mod,mod1)
anova(mod,mod1)
# Major.Lineage appears necessary (but oddly enough with 1-yr bandwidth, this suggests it's ok to leave out the Major.Lineage data)

# let's examine cluster with economic indicators and Major.lineage:
mod_clust1 <- glm(seqden~fridge+carmotobike+bed+radio+phone+Major.Lineage+cluster,dat,family=Gamma(link='identity'))
mod_clust0 <- glm(seqden~fridge+carmotobike+bed+radio+phone+Major.Lineage,dat[!is.na(dat$cluster),],family=Gamma(link='identity'))
AIC(mod_clust0,mod_clust1)

#not favored according to AIC... could it be a proxy for economic status? 
mod_clust1 <- glm(seqden~Major.Lineage+cluster,dat,family=Gamma(link='identity'))
mod_clust0 <- glm(seqden~Major.Lineage,dat[!is.na(dat$cluster),],family=Gamma(link='identity'))
AIC(mod_clust0,mod_clust1)

# just a marginal improvement over Major.Lineage, ~8 AIC values.. so borderline


# The missingness patterns are kind of tricky.
# It would be nice to include the other economic variables that got tossed out for containing missingness. It seems 
# like they all have 57 missing values - not too bad

# clinic/hospital visit questions seem to all have 113 missing entries

# Others are much wosre - last_worked_hs has 695 missing entries



# lump the economic indicators together - haselectricity has 113 NA, the others 57 or 59, if any
economic <- c('fridge','carmotobike','bed','radio','phone','peopleinhh','sleepinsameroom','cookinglocation','smoke','knowanyonewithtb','haselectricity')
# fridge, carmotobike, bed, radio, and phone are all positively correlated look at fisher.test(with(dat,table(fridge,carmotobike)), etc.)


# fit the model to the smaller set of data that has information about all of the economic indicators:
mod2 <- glm(seqden~fridge+carmotobike+bed+radio+phone+peopleinhh+sleepinsameroom+cookinglocation+smoke+knowanyonewithtb+haselectricity+Major.Lineage+culture_positive,dat,family=Gamma(link='identity'))

# do the additional economic indicators carry any extra information than the old ones did?
crud <- apply(dat[,economic],1,is.na)
crud <- t(crud)
crud <- apply(crud, 1, sum)
mod2o <- glm(seqden~fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive,dat[crud==0,],family=Gamma(link='identity'))
rm(crud)
AIC(mod2,mod2o)
anova(mod2,mod2o)
# doesn't look like the extra economic information is telling us much (but we are fitting to fewer data points now)

# are these newer economic indicators just telling us what we already knew?
mod2a <- glm(seqden~peopleinhh+sleepinsameroom+cookinglocation+smoke+knowanyonewithtb+haselectricity+Major.Lineage+culture_positive,dat,family=Gamma(link='identity')
)
AIC(mod2,mod2a)
anova(mod2,mod2a)
# actually, it looks like they're less useful than the older economic data (fridge, carmotobike, bed, radio, phone)



# let's see if anyhivclinic and anyop add anything:
modhc <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler+anyhivclinic+anyop,dat,family=Gamma(link='identity'))
AIC(modhc,mod)
anova(modhc,mod)
# doesn't look like it, but maybe the more detailed information says something interesting...

# let's see if lastclinicvisit does anything:
mod_clin0 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler+anyhivclinic+anyop,dat[!is.na(dat$lastclinicvisit),],family=Gamma(link='identity'))

mod_clin1 <- glm(seqden~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler+anyhivclinic+anyop+lastclinicvisit,dat,family=Gamma(link='identity'))

#convergence issues... not sure why

# the vast majority of sequences had their last clinic visit in the last 6 months.
# there isn't an obvious association between LBI and lastclinicvisit. LBI might be skewed high...
# but it seems difficult to tell

table(dat$lastclinicvisit)
boxplot(seqden~lastclinicvisit,dat)

# random forest
require(randomForest)


# first let's set a threshold for high density
dat$highseqden <- 0
dat$highseqden[dat$seqden > 5] <- 1

#randomForest will break if there are columns in the data with all NA's; just select the columns we want

rfdat <- dat[,c('x04fac_code','x05year','x09iptpst','x10iptdiag','x17tbtype','x16patcat','x20hiv','x21art','x22cotri','x02res','sex','age','smeartest','sympcough','sympsweat','sympfever','sympweight','sympblood','sympbreath','fridge','carmotobike','bed','radio','phone','Major.Lineage','Drug.resistance.Tbprofiler','anyhivclinic','anyop','seqden')] 

#split rf dat into a training and test set
set.seed(117)
ind <- sample(2, nrow(rfdat), replace=T, prob=c(0.7, 0.3))
train <- rfdat[ind==1,]
test <- rfdat[ind==2,]

rf <- randomForest(seqden~., rfdat) 

par(mfrow=c(1,2))
p1 <- predict(rf, train)
plot(train$seqden,p1, xlab='LBI (training set)', ylab='Predicted LBI', main='Training (70%)',
	xlim=c(0,20),ylim=c(0,20))
abline(a=0,b=1)

p2 <- predict(rf, test)
plot(test$seqden, p2, xlab='LBI (test set)', ylab='Predicted LBI', main='Test (30%)',
	xlim=c(0,20), ylim=c(0,20))
abline(a=0,b=1)

varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf)

# not totally out of line with the glm results, although it does think x04fac_code and age are quite important

# 551 datapoints have a value for cluster (a nbd identifier)
# let's see how including this improves things:

rfdat_clust <- dat[!is.na(dat$clust),c('x04fac_code','x05year','x09iptpst','x10iptdiag','x17tbtype','x16patcat','x20hiv','x21art','x22cotri','x02res','sex','age','smeartest','sympcough','sympsweat','sympfever','sympweight','sympblood','sympbreath','fridge','carmotobike','bed','radio','phone','Major.Lineage','Drug.resistance.Tbprofiler','anyhivclinic','anyop','seqden','cluster')] 

#split rf dat into a training and test set
set.seed(117)
ind <- sample(2, nrow(rfdat_clust), replace=T, prob=c(0.7, 0.3))
train <- rfdat_clust[ind==1,]
test <- rfdat_clust[ind==2,]

rf_clust <- randomForest(seqden~., rfdat_clust) 

par(mfrow=c(1,2))
p1 <- predict(rf_clust, train)
plot(train$seqden,p1, xlab='LBI (training set)', ylab='Predicted LBI', main='Training (70%)',
	xlim=c(0,20),ylim=c(0,20))
abline(a=0,b=1)

p2 <- predict(rf_clust, test)
plot(test$seqden, p2, xlab='LBI (test set)', ylab='Predicted LBI', main='Test (30%)',
	xlim=c(0,20), ylim=c(0,20))
abline(a=0,b=1)

par(mfrow=c(1,1))
varImpPlot(rf_clust,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_clust)

# cluster is identified as very important
boxplot(seqden~cluster,dat,xlab='cluster',ylab='LBI')

# let's try including lastclinicvisit in random forest:
rfdat_lastclinicvisit <- dat[!is.na(dat$lastclinicvisit),c('x04fac_code','x05year','x09iptpst','x10iptdiag','x17tbtype','x16patcat','x20hiv','x21art','x22cotri','x02res','sex','age','smeartest','sympcough','sympsweat','sympfever','sympweight','sympblood','sympbreath','fridge','carmotobike','bed','radio','phone','Major.Lineage','Drug.resistance.Tbprofiler','anyhivclinic','anyop','seqden','lastclinicvisit')] 

#split rf dat into a training and test set
set.seed(117)
ind <- sample(2, nrow(rfdat_lastclinicvisit), replace=T, prob=c(0.7, 0.3))
train <- rfdat_lastclinicvisit[ind==1,]
test <- rfdat_lastclinicvisit[ind==2,]

rf_lastclinicvisit <- randomForest(seqden~., rfdat_lastclinicvisit) 

par(mfrow=c(1,2))
p1 <- predict(rf_lastclinicvisit, train)
plot(train$seqden,p1, xlab='LBI (training set)', ylab='Predicted LBI', main='Training (70%)',
	xlim=c(0,20),ylim=c(0,20))
abline(a=0,b=1)

p2 <- predict(rf_lastclinicvisit, test)
plot(test$seqden, p2, xlab='LBI (test set)', ylab='Predicted LBI', main='Test (30%)',
	xlim=c(0,20), ylim=c(0,20))
abline(a=0,b=1)

par(mfrow=c(1,1))
varImpPlot(rf_lastclinicvisit,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_lastclinicvisit)

# lastclinicvisit doesn't make the top ten 
boxplot(seqden~lastclinicvisit,dat,xlab='lastclinicvisit',ylab='LBI')



# try a lasso regression on those same variables
require(glmnet)
lassodat <- rfdat

x_vars <- model.matrix(seqden~. , lassodat)
y_var <- lassodat$seqden
lambda_seq <- 10^seq(100, -2, by= -.1)

#split test data into test and train
set.seed(117)
train = sample(1:nrow(x_vars), nrow(x_vars)/2)
x_test = (-train)

cv_output <- cv.glmnet(x_vars[train,], y_var[train],
		alpha = 1, lambda = lambda_seq,
		nfolds = 5)

best_lam <- cv_output$lambda.min
best_lam # blows up

lambda <- .2
lasso_best <- glmnet(x_vars[train,], y_var[train],
		alpha = 1, lambda = lambda)
coef(lasso_best)

# penalized Maximum Likelihood wants to select lambda=infinty apparently (?), but
# coefficients get shrunk to zero roughly according to their level importance we 
# saw using AIC  

# kind of weird that lasso wants to kill coefficients for particular levels of a factor. 
# Let's try a group lasso to get rid of all levels of a factor at once.
require(gglasso)

# introduce a design matrix

modmat <- model.matrix(~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler,dat)

# now, define groups for the group lasso:
groups <- c(
	1,	
	rep(2,11),
	3,
	4,
	5,
	6,
	7,7,7,7,
	8,8,
	9,
	10,
	11,
	12,
	13,
	14,
	15,15,15,15,15,15,
	16,16,16,16,16,
	17,17,17,
	18,
	19,19)

gr <- gglasso(modmat,dat$seqden,lambda=exp(seq(-1,-10,length=30)), group=groups, loss='ls',intercept=F)
#coef(gr)

# plot the coefficient trajectories (with labels):

beta <- gr$beta
lambdas <- gr$lambda

plot_data <- data.frame(
	lambda=rep(log(lambdas),each=nrow(beta)),
	coefficient = as.vector(beta),
	variable = rep(rownames(beta),length(lambdas))
)

lasplot <- ggplot(plot_data, aes(x = lambda, y = coefficient, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Coefficients vs. log(Lambda)",
       x = "log(Lambda)", y = "Coefficients") +
  geom_text(
    data = plot_data[plot_data$lambda == min(log(lambdas)), ], 
    aes(label = variable), 
    hjust = 1.1, size = 3
  ) +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "lines"))+  # Increase right margin
  xlim(-12,0)

lasplot

#load in the BLAST.xls file with "cluster", a nbd identifier
require(readxl)
blast <- read_excel('BLAST.xls')

# each of dat and blast has pid's that are not in the other (551 in common)
dat <- merge(dat, blast[,c('pid','cluster')], all.x=T)

# let's rerun the lasso to see if it thinks cluster is important

# introduce a design matrix
modmat <- model.matrix(~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler+cluster,dat)

# now, define groups for the group lasso:
groups <- c(
	1,	
	rep(2,11),
	3,
	4,
	5,
	6,
	7,7,7,
	8,8,
	9,
	10,
	11,
	12,
	13,
	14,
	15,15,15,15,15,15,
	16,16,16,16,16,
	17,17,17,
	18,
	19,19,
	rep(20,72))

gr <- gglasso(modmat,dat[!is.na(dat$cluster),]$seqden,lambda=exp(seq(-1,-10,length=30)), group=groups, loss='ls',intercept=F)
#coef(gr)

# plot the coefficient trajectories (with labels):

beta <- gr$beta
lambdas <- gr$lambda

plot_data <- data.frame(
	lambda=rep(log(lambdas),each=nrow(beta)),
	coefficient = as.vector(beta),
	variable = rep(rownames(beta),length(lambdas))
)

lasplot <- ggplot(plot_data, aes(x = lambda, y = coefficient, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Coefficients vs. log(Lambda)",
       x = "log(Lambda)", y = "Coefficients") +
  geom_text(
    data = plot_data[plot_data$lambda == min(log(lambdas)), ], 
    aes(label = variable), 
    hjust = 1.1, size = 3
  ) +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), "lines"))+  # Increase right margin
  xlim(-12,0)

lasplot

# now, we want to look at healthcare utilization
# lastclinicvisit records the last time they visited a clinic or hospital for any reason
# a small number of people answered questions about the timing of their last 3 hospital admissions
# 	(admissionunittone (n=123), admissionunittwo (n=53), admissionunitthree (n=40))

# let's just see if we can see anything from lastclinicvisit

# introduce a design matrix
modmat <- model.matrix(~x04fac_code+x05year+x09iptpst+x10iptdiag+x17tbtype+x16patcat+x20hiv+x21art+x22cotri+x02res+sex+age+smeartest+sympcough+sympsweat+sympfever+sympweight+sympblood+sympbreath+fridge+carmotobike+bed+radio+phone+Major.Lineage+culture_positive+Drug.resistance.Tbprofiler+cluster+lastclinicvisit,dat)

# now, define groups for the group lasso:
groups <- c(
	1,	
	rep(2,10), # one of the locations is not represented in this subset of the data
	3,
	4,
	5,
	6,
	7,7,7,
	8,8,
	9,
	10,
	11,
	12,
	13,
	14,
	15,15,15,15,15,15,
	16,16,16,16,16,
	17,17,17,
	18,
	19,19,
	rep(20,72),
	rep(21,4))

#gr <- gglasso(modmat,dat[!is.na(dat$cluster) & !is.na(dat$lastclinicvisit),]$seqden,lambda=exp(seq(-1,-10,length=30)), group=groups, loss='ls',intercept=F)


# plot the coefficient trajectories (with labels):

beta <- gr$beta
lambdas <- gr$lambda

plot_data <- data.frame(
        lambda=rep(log(lambdas),each=nrow(beta)),
        coefficient = as.vector(beta),
        variable = rep(rownames(beta),length(lambdas))
)

lasplot <- ggplot(plot_data, aes(x = lambda, y = coefficient, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Coefficients vs. log(Lambda)",
       x = "log(Lambda)", y = "Coefficients") +
  geom_text(
    data = plot_data[plot_data$lambda == min(log(lambdas)), ],
    aes(label = variable),
    hjust = 1.1, size = 3
  ) +
  theme(legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "lines"))+  # Increase right margin
  xlim(-12,0)

lasplot


# some ideas
# 1. think about a structured coalescent model with 2 hidden populations; can we detect the presence of pop structure?
# 2. testing might be biased toward particular communities. It would be nice to recover differential prevalence or incidence 
# 3. Or maybe ask: Given that we think our case-finding is targeted towards a particular group, we would like to use
#    the genomic data to help us understand how much of the transmission is actually occurring within that group vs. 
#    migrating into it

#
d4 <- cophenetic(timetree4)

# if we stratify by hivstatus, do sequences look more or less similar within the positive and negative groups?
crud = dat[dat$Major.Lineage=='lineage4',]

crud = crud[,c('Sequence_name','hivstatus')]
hivpos = crud[crud$hivstatus=='Positive',1]
hivneg = crud[crud$hivstatus=='Negative',1]

# write the numeric indices for pos/neg
inds.hivpos <- unlist(sapply(hivpos, function(x) which(rownames(d4)==x)))
inds.hivpos <- sort(unname(inds.hivpos))

inds.hivneg <- unlist(sapply(hivneg, function(x) which(rownames(d4)==x)))
inds.hivneg <- sort(unname(inds.hivneg))


d4.hivpos <- d4[inds.hivpos,inds.hivpos]
d4.hivneg <- d4[inds.hivneg,inds.hivneg]

d4.hivstatus <- d4[c(inds.hivpos,inds.hivneg),c(inds.hivpos,inds.hivneg)]

# within-hivpositive distances:
dij.hivpos <- d4.hivpos
dij.hivpos <- dij.hivpos[upper.tri(dij.hivpos)]

# within-hivnegitive distances:
dij.hivneg <- d4.hivneg
dij.hivneg <- dij.hivneg[upper.tri(dij.hivneg)]

# across-group distances:
dij.cross <- d4[inds.hivpos,inds.hivneg]

hivneg <- data.frame(hivstatus='hivneg',dist=dij.hivneg)
hivpos <- data.frame(hivstatus='hivpos',dist=dij.hivpos)
hivcross <- data.frame(hivstatus='hivcross',dist=c(dij.cross))

hivdat <- rbind(hivneg,hivpos,hivcross)

# doesn't look like distances within or between hivstatus groups are different:
plot(dist~as.factor(hivstatus),hivdat)

#what about "phone"?

crud = dat[dat$Major.Lineage=='lineage4',]

crud = crud[,c('Sequence_name','phone')]
phoneyes = crud[crud$phone=='Yes',1]
phoneno = crud[crud$phone=='No',1]

# write the numeric indices for pos/neg
inds.phoneyes <- unlist(sapply(phoneyes, function(x) which(rownames(d4)==x)))
inds.phoneyes <- sort(unname(inds.phoneyes))

inds.phoneno <- unlist(sapply(phoneno, function(x) which(rownames(d4)==x)))
inds.phoneno <- sort(unname(inds.phoneno))


d4.phoneyes <- d4[inds.phoneyes,inds.phoneyes]
d4.phoneno <- d4[inds.phoneno,inds.phoneno]

d4.hivstatus <- d4[c(inds.phoneyes,inds.phoneno),c(inds.phoneyes,inds.phoneno)]

# within-phoneyesitive distances:
dij.phoneyes <- d4.phoneyes
dij.phoneyes <- dij.phoneyes[upper.tri(dij.phoneyes)]

# within-phonenoitive distances:
dij.phoneno <- d4.phoneno
dij.phoneno <- dij.phoneno[upper.tri(dij.phoneno)]

# across-group distances:
dij.cross <- d4[inds.phoneyes,inds.phoneno]

phoneno <- data.frame(phone='phoneno',dist=dij.phoneno)
phoneyes <- data.frame(phone='phoneyes',dist=dij.phoneyes)
phonecross <- data.frame(phone='phonecross',dist=c(dij.cross))

phonedat <- rbind(phoneno,phoneyes,phonecross)

# doesn't look like distances within or between phone groups are different:
plot(dist~as.factor(phone),phonedat)


