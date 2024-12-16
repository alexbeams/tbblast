rm(list=ls())

# Calculate delta-statistics based on the code here:
# https://github.com/mrborges23/delta_statistic/blob/master/README.md

library(parallel)
library(ape)
source('delta_statistic/code.R')

# let's not worry about estimating the parameters of the entropy distribution yet - instead,
# we calculate state entropies for each of the categorical variables in the data
# randomize the traits to calculate entropies for comparison
 
#load in the lineage 1-3 posteriors from BEAST:
lin1phi <- read.nexus("lineage_files/beast_files/alignment1.beauti-alignment_lineage1.trees")
lin2phi <- read.nexus("lineage_files/beast_files/alignment2.beauti-alignment_lineage2.trees")
lin3phi <- read.nexus("lineage_files/beast_files/alignment3.beauti-alignment_lineage3.trees")

# load in the lineage 4 posterior from delphy:



# load in the metadata and the lineage/drug resistance data
drugdat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')
names(drugdat)[1] <- 'Sequence_name'

# merge them by sequence name
dat <- merge(blastdat, drugdat, by='Sequence_name')

# load in the dataframe that contains "cluster", a nbd identifier
require(readxl)
blast <- read_excel('BLAST.xls')

# each of dat and blast has pid's that are not in the other (551 in common)
dat <- merge(dat, blast[,c('pid','cluster')], all.x=T)

dat1 <- dat[dat$Major.Lineage=='lineage1',]
dat2 <- dat[dat$Major.Lineage=='lineage2',]
dat3 <- dat[dat$Major.Lineage=='lineage3',]
dat4 <- dat[dat$Major.Lineage=='lineage4',]


#subsample from the posterior samples:
nsamps <- 100
phi1_sample <- sample(lin1phi,size=nsamps)
phi2_sample <- sample(lin2phi,size=nsamps)
phi3_sample <- sample(lin3phi,size=nsamps)



## try calculating delta for the first tree
#
#tree <- phi1_sample[[1]]
#
## need to ensure that all branches have positive length; if any are zero, just jitter them a bit:
#tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
#
## try 'phone' as the trait to start:
#trait <- dat[dat$Major.Lineage=='lineage1','hivstatus']
#
#trait <- as.numeric(trait=='Positive')
#trait[is.na(trait)] <- 0
#names(trait) <- dat[dat$Major.Lineage=='lineage1','Sequence_name']
#
# this is not that slow, and we can calculate it in parallel across the posterior trees:
#deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100)

# this is slow
#random_delta <- rep(NA,100)
#for (i in 1:100){
#  rtrait <- sample(trait)
#  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
#}
#p_value <- sum(random_delta>deltaA)/length(random_delta)
#boxplot(random_delta)
#abline(h=deltaA,col="red")

getDeltas <- function(posterior_tree_sample,trait) mclapply(posterior_tree_sample, function(x) delta(trait,x,0.1,.06,10000,10,100),mc.cores=8) 

getTrait <- function(traitnm,lineage){
	trait = dat[dat$Major.Lineage==lineage,traitnm]
	names(trait) <- dat[dat$Major.Lineage==lineage,'Sequence_name']	
	trait <- factor(trait)
	return(trait)
}

getEntropy <- function(posterior_tree, trait,randomize=F){
	if(randomize==T){crud=trait[sample(1:length(trait))];names(crud)<-names(trait);trait<-crud}	
	ans <- ace(trait,posterior_tree,type='d')
	ar <- nentropy(ans$lik.anc)
	return(mean(ar))

}

getEntropies <- function(posterior_tree_sample, trait,randomize=F) mclapply(posterior_tree_sample, function(x) getEntropy(x,trait,randomize=randomize) , mc.cores=8)

plotAncRes <- function(posterior_tree,trait,co){
	#trait <- getTrait(traitnm,lineage)	
	plot(posterior_tree, type='fan', FALSE, label.offset=0, cex=0.1)
	tiplabels(pch=22, bg=co[as.numeric(trait)], cex=1)
	ans <- ace(trait,posterior_tree,type='d')
	nodelabels(thermo=ans$lik.anc,piecol=co,cex=0.3)	

}

varnms <- c('hivstatus','smoke','sex','radio','phone','fridge','sympcough','tbcategory','tbclass')

x <- lapply(varnms, function(x) getEntropies(phi1_sample,getTrait(x,'lineage1'))  )

par(mfrow=c(3,3))
for(i in 1:9){ hist(unlist(x[[i]]),xlim=c(0,1),main=varnms[i])  }

# most of the entropies seem distributed close to 0 or 1
# fridge seems interesting because it seems further from the boundaries


plotEntropies <- function(posterior_sample,traitnm,linnm){
	x <- getEntropies(posterior_sample,getTrait(traitnm,linnm))
	y <- getEntropies(posterior_sample,getTrait(traitnm,linnm),randomize=T)
	xmin = 0.5*min(c(unlist(x),unlist(y)))	
	xmax = min(1,1.5*max(c(unlist(x),unlist(y)))	)

	hist(unlist(y),col='red',xlim=c(xmin,xmax),main=paste(traitnm,linnm))
	hist(unlist(x),add=T,col='blue')
}

par(mfrow=c(1,1))
plotAncRes(phi1_sample[[1]],getTrait('fridge','lineage1'),c('red','blue'))

#compare the counts of fridge and phone:
table(dat[dat$Major.Lineage=='lineage1','fridge'])
table(dat[dat$Major.Lineage=='lineage1','phone'])

#what happens if we randomly change the fridge variable?

trait = getTrait('fridge','lineage1')
crud = sample(trait)
names(crud) <- names(trait)

fridge_rand <- getEntropies(phi1_sample,crud)

# looks almost identical to the entropy distribution for fridge; numbers must matter a lot


