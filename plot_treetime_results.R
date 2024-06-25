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

#generate local branching indices for the tips
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
dat1.merge <- merge(data.frame(Sequence_name=timetree1$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat2.merge <- merge(data.frame(Sequence_name=timetree2$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat3.merge <- merge(data.frame(Sequence_name=timetree3$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

dat4.merge <- merge(data.frame(Sequence_name=timetree4$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

datMbovis.merge <- merge(data.frame(Sequence_name=timetreeMbovis$tip.label), dat,
	by.x="Sequence_name",by.y="Sequence_name",all.x=TRUE )

# lineage-specific ggtree objects:
p1.dat <- ggtree(timetree1,layout='circular')
p2.dat <- ggtree(timetree2,layout='circular')
p3.dat <- ggtree(timetree3,layout='circular')
p4.dat <- ggtree(timetree4,layout='circular')
pMbovis.dat <- ggtree(timetreeMbovis,layout='circular')



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



# are any metadata variables correlated with LBI?

#any noticeable trends in time?
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

