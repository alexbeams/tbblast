# Clear the environment
rm(list=ls())

# one can choose bandwidth parameters for LBI, THD, and inverse-exponential weighted average path lengths
#	to make them all basically the same (correlations close to 80%)

# How do these values compare to the likelihoods of the individual sequences? Would that really be telling us
#	anything about fitness?

# we care about high levels of relatedness as a proxy for recent transmission/higher fitness


# Load necessary libraries
library(ape) #tree stuff
library(ggtree)
library(ggplot2)
library(gridExtra) #plot multiple ggtrees in windows
library(igraph) #calculate centrality values


source('lbi.R') #load in bdearlove's lbi function

treelayout <- 'circular'
dotsize <- 1 

#set.seed(4927)

# Generate a random tree with numtips tips
numtips <- 600
tree <- rtree(numtips)
tree$tip.label <- paste0('t', 1:length(tree$tip.label))
tree$node.label <- paste0('n', 1:tree$Nnode)

# Calculate path lengths between all tips
tip_lengths <- cophenetic(tree)

# Calculate path lengths between all of the internal nodes, as well as the tips
# numtips x numtips upper-left block of the dist.nodes() output equals the cophenetic output
allnode_lengths <- dist.nodes(tree)

# Compute the index as the average path length to other tips only
tip.mean <- apply(tip_lengths,1,mean) 
tip.var <- apply(tip_lengths,1,var)
tip.precision <- 1/apply(tip_lengths,1,var)
tip.min <- apply(tip_lengths,1, function(x) min(x[x>0]))
tip.recipmax <- 1/apply(tip_lengths,1,max)
tip.max <- apply(tip_lengths,1,max)
tip.prox <- apply(tip_lengths,1,function(x) mean(1/x[x>0]) ) 
#  mean patristic distance overall to parameterize lbi and exp weighting scheme:
tau <- mean(tip_lengths[upper.tri(tip_lengths)],diag=F)
tip.exp <- apply(tip_lengths, 1, function(x) sum(exp(-x/.0625/tau))  ) 



# Create a data frame with the mean tip-to-tip lengths 
df.tip <- data.frame(tip = tree$tip.label, 
	tip.mean = tip.mean,
	tip.var = tip.var,
	tip.precision = tip.precision,
	tip.min = tip.min,
	tip.recipmax = tip.recipmax,
	tip.max = tip.max,
	tip.prox = tip.prox,
	tip.exp = tip.exp,
	stringsAsFactors = FALSE)

#calculate the THD using the tip distances first:
source('thd.R')

m <- 1500 #46000
mu <- 5e-4
timescale <- .2 
bandwidth <- tmrca2bandwidth(timescale, m, mu) 
df.tip$tip.thd <- thd(tip_lengths, bandwidth, m)

# Compute the as the average patristic distance to internal nodes as well as tips
allnode.mean <- apply(allnode_lengths,1,mean) 
allnode.var <- apply(allnode_lengths,1,var)
allnode.precision <- 1/allnode.var 
allnode.min <- apply(allnode_lengths,1, function(x) min(x[x>0]))
allnode.recipmax <- 1/apply(allnode_lengths,1,max)
allnode.max <- apply(allnode_lengths,1,max)
allnode.prox <- apply(allnode_lengths, 1, function(x) mean(1/x[x>0.01])  ) 
#  mean patristic distance overall to parameterize lbi and exp weighting scheme:
tau <- mean(allnode_lengths[upper.tri(allnode_lengths)],diag=F)
allnode.exp <- apply(allnode_lengths, 1, function(x) sum(exp(-x/.0625/tau))  ) 
allnode.lbi <- lbi(tree,tau=0.0625 * tau)


# Create a ggtree plot with a circular layout
p <- ggtree(tree, layout = treelayout)

# Add the data frame to the plot
p.tip <- p %<+% df.tip
 
p.tip.mean <- p.tip + geom_tippoint(aes(color = tip.mean),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Mean Patristic Distance to Other Tips')

p.tip.var <- p.tip + geom_tippoint(aes(color = tip.var),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Variance in Patristic Distance to other tips')

p.tip.precision <- p.tip + geom_tippoint(aes(color = tip.precision),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Precision in Patristic Distance to other tips')

p.tip.min <- p.tip + geom_tippoint(aes(color = tip.min),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Minimum Patristic Distance to other tips')

p.tip.max <- p.tip + geom_tippoint(aes(color = tip.max),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Max Patristic Distance to other tips')

p.tip.recipmax <- p.tip + geom_tippoint(aes(color = tip.recipmax),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Reciprocal Max Patristic Distance to other tips')

p.tip.exp <- p.tip + geom_tippoint(aes(color = tip.exp),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Exp weighted Patristic Distance to other tips')

p.tip.prox <- p.tip + geom_tippoint(aes(color = tip.prox),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Proximity (1/distance) between tips')

p.tip.thd <- p.tip + geom_tippoint(aes(color = tip.thd),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('THD relative to tips only')


# Create a data frame with the mean allnode-to-allnode lengths 
df.allnode <- data.frame(allnode = c(tree$tip.label,tree$node.label), 
	allnode.mean = allnode.mean, 
	allnode.var = allnode.var,
	allnode.precision = allnode.precision,	
	allnode.min = allnode.min,
	allnode.exp = allnode.exp,
	allnode.lbi = allnode.lbi,
	allnode.max = allnode.max,
	allnode.recipmax = allnode.recipmax,
	allnode.prox = allnode.prox,
	stringsAsFactors = FALSE)

# try calculating THD using the patristic distances for all the nodes
df.allnode$allnode.thd <- thd(allnode_lengths, bandwidth, m)

# convert the tree to a graph, to calculate some centrality values
edges <- tree$edge

# Get the node labels (tip labels and internal node labels)
node_labels <- c(tree$tip.label, tree$node.label)

# Create the igraph object
g <- graph_from_edgelist(as.matrix(edges), directed = FALSE)

# Add node labels as vertex attributes
V(g)$name <- node_labels

# Add the lengths as edge attributes
E(g)$weight <- tree$edge.length 

allnode.betweenness1 <- betweenness(g,cutoff=1,)
allnode.betweenness2 <- betweenness(g,cutoff=2,)
allnode.betweenness4 <- betweenness(g,cutoff=4,)
allnode.betweenness8 <- betweenness(g,cutoff=8,)
allnode.betweenness16 <- betweenness(g,cutoff=16,)

df.allnode$allnode.betweenness1 <- allnode.betweenness1
df.allnode$allnode.betweenness2 <- allnode.betweenness2
# ^ this one basically seems to be measuring clade size
df.allnode$allnode.betweenness4 <- allnode.betweenness4
df.allnode$allnode.betweenness8 <- allnode.betweenness8
df.allnode$allnode.betweenness16 <- allnode.betweenness16

#some others:
df.allnode$allnode.pagerank <- page_rank(g)[[1]]
df.allnode$allnode.diversity <- diversity(g)
df.allnode$allnode.alphacentrality <- alpha_centrality(g)
df.allnode$allnode.eigencentrality <- eigen_centrality(g)[[1]]
df.allnode$allnode.powercentrality <- power_centrality(g)



# Create a ggtree plot with a circular layout
p <- ggtree(tree, layout = treelayout)

# Add the data frame to the plot
p.allnode <- p %<+% df.allnode
 
p.allnode.mean <- p.allnode + geom_tippoint(aes(color = allnode.mean),size=dotsize) +
	geom_nodepoint(aes(color = allnode.mean),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Mean Patristic Distance to all other nodes')


#grid.arrange(p.tip, p.allnode.mean, ncol = 2)

# These look good - basically the outside of the trees are "warmest", and using mean distance between tips
#  is basically the same as mean distance between all nodes - but maybe this could change if we were using
#  some other, more imbalanced way of producing trees

p.allnode.var <- p.allnode + geom_tippoint(aes(color = allnode.var),size=dotsize) +
	geom_nodepoint(aes(color = allnode.var),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Variance in Patristic Distance to all other nodes')

p.allnode.precision <- p.allnode + geom_tippoint(aes(color = allnode.precision),size=dotsize) +
	geom_nodepoint(aes(color = allnode.precision),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Precision in Patristic Distance to all other nodes')

p.allnode.min <- p.allnode + geom_tippoint(aes(color = allnode.min),size=dotsize) +
	geom_nodepoint(aes(color = allnode.min),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Minimum Patristic Distance to all other nodes')

p.allnode.max <- p.allnode + geom_tippoint(aes(color = allnode.max),size=dotsize) +
	geom_nodepoint(aes(color = allnode.max),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Max Patristic Distance to all other nodes')

p.allnode.recipmax <- p.allnode + geom_tippoint(aes(color = allnode.recipmax),size=dotsize) +
	geom_nodepoint(aes(color = allnode.recipmax),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Reciprocal Max Patristic Distance to all other nodes')

p.allnode.exp <- p.allnode + geom_tippoint(aes(color = allnode.exp),size=dotsize) +
	geom_nodepoint(aes(color = allnode.exp),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Exp weighted Patristic Distance to all other nodes')

p.allnode.lbi <- p.allnode + geom_tippoint(aes(color = allnode.lbi),size=dotsize) +
	geom_nodepoint(aes(color = allnode.lbi),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('LBI (requires all nodes)')

p.allnode.prox <- p.allnode + geom_tippoint(aes(color = allnode.prox),size=dotsize) +
	geom_nodepoint(aes(color = allnode.prox),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Proximity (1/distance) to all other nodes')

p.allnode.thd <- p.allnode + geom_tippoint(aes(color = allnode.thd),size=dotsize) +
	geom_nodepoint(aes(color = allnode.thd),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('THD relative to all other nodes')

p.allnode.betweenness1 <- p.allnode + #geom_tippoint(aes(color = allnode.betweenness1),size=dotsize) +
	geom_nodepoint(aes(color = allnode.betweenness1),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Betweenness centrality (cutoff of 1 patristic distance units)')

p.allnode.betweenness2 <- p.allnode + #geom_tippoint(aes(color = allnode.betweenness2),size=dotsize) +
	geom_nodepoint(aes(color = allnode.betweenness2),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Betweenness centrality (cutoff of 2 patristic distance units)')

p.allnode.betweenness4 <- p.allnode + #geom_tippoint(aes(color = allnode.betweenness4),size=dotsize) +
	geom_nodepoint(aes(color = allnode.betweenness4),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Betweenness centrality (cutoff of 4 patristic distance units)')

p.allnode.betweenness8 <- p.allnode + #geom_tippoint(aes(color = allnode.betweenness8),size=dotsize) +
	geom_nodepoint(aes(color = allnode.betweenness8),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Betweenness centrality (cutoff of 8 patristic distance units)')

p.allnode.betweenness16 <- p.allnode + #geom_tippoint(aes(color = allnode.betweenness16),size=dotsize) +
	geom_nodepoint(aes(color = allnode.betweenness16),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Betweenness centrality (cutoff of 16 patristic distance units)')

p.allnode.pagerank <- p.allnode + #geom_tippoint(aes(color = allnode.pagerank),size=dotsize) +
	geom_nodepoint(aes(color = allnode.pagerank),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('PageRank')


p.allnode.eigencentrality <- p.allnode + #geom_tippoint(aes(color = allnode.eigencentrality),size=dotsize) +
	geom_nodepoint(aes(color = allnode.eigencentrality),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Eigen Centrality')

p.allnode.alphacentrality <- p.allnode + #geom_tippoint(aes(color = allnode.alphacentrality),size=dotsize) +
	geom_nodepoint(aes(color = allnode.alphacentrality),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Alpha Centrality')

p.allnode.powercentrality <- p.allnode + #geom_tippoint(aes(color = allnode.powercentrality),size=dotsize) +
	geom_nodepoint(aes(color = allnode.powercentrality),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Power Centrality')

p.allnode.diversity <- p.allnode + #geom_tippoint(aes(color = allnode.diversity),size=dotsize) +
	geom_nodepoint(aes(color = allnode.diversity),size=dotsize) +
        scale_color_gradient(high = 'red', low = 'blue') +
        theme_tree2() +
        ggtitle('Diversity' )


# this looks really similar to p.allnode.exp
#print(p.allnode.lbi)

#hist(lbi(tree,tau=.0625*tau)) 
#the value of .0625*mean(patristic distance) gives 
# a distribution of lbi that looks pretty close to an exponential, perhaps?

#print(p.tip.prox)
grid.arrange(p.allnode.exp, p.allnode.lbi, p.allnode.thd,
	p.allnode.betweenness1, ncol = 2,nrow=2)

cor1 <- cor(df.allnode$allnode.exp,df.allnode$allnode.lbi)
cor2 <- cor(df.allnode$allnode.exp,df.allnode$allnode.thd)
cor3 <- cor(df.allnode$allnode.lbi,df.allnode$allnode.thd)
cortab <- c(cor1,cor2,cor3)
names(cortab) <- c('Exp-LBI','Exp-THD','LBI-THD')

print('correlations: ')
print(cortab)

grid.arrange(p.allnode.betweenness1, p.allnode.betweenness2, p.allnode.betweenness4,
	p.allnode.betweenness8, ncol = 2,nrow=2)

