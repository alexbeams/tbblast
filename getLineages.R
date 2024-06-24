rm(list=ls())

# Identify and extract lineages from the IQ-TREE output. We will generate
# timed trees on the subtrees consisting of the individual lineages

require(ape)
require(phytools)


# read in the .nexus timetree generated from the iq-tree kimura model
# as well as the data file from Ben S. with lineage info

# read in the un-timed tree:
tree <- read.tree('constant_site_correction_kimura/Malawi_final_filtered.treefile')

lineagedat <- read.csv('Malawi_final_stats.csv')
blastdat <- read.csv('BLAST_ePAL_SeqID_NoGPS.csv')

#there are some lineages in lineagedat that are not included in blastdat, 
# i.e. not included in tree

extras <- setdiff(lineagedat$Sequence.name, blastdat$Sequence_name)
extrainds <- c()
for(i in 1:length(extras)) extrainds[i] <- which(lineagedat$Sequence.name == extras[i]) 

lineagedat <- lineagedat[-extrainds,]

# Identify the lineages:

lineage1_tips <- lineagedat[c(lineagedat$Major.Lineage=='lineage1'),'Sequence.name']
lineage2_tips <- lineagedat[c(lineagedat$Major.Lineage=='lineage2'),'Sequence.name']
lineage3_tips <- lineagedat[c(lineagedat$Major.Lineage=='lineage3'),'Sequence.name']
lineage4_tips <- lineagedat[c(lineagedat$Major.Lineage=='lineage4'),'Sequence.name']
lineageMbovis_tips <- lineagedat[c(lineagedat$Major.Lineage=='M.bovis'),'Sequence.name']

# Extract subtrees
lineage1_tree <- keep.tip(tree, lineage1_tips)
lineage2_tree <- keep.tip(tree, lineage2_tips)
lineage3_tree <- keep.tip(tree, lineage3_tips)
lineage4_tree <- keep.tip(tree, lineage4_tips)
lineageMbovis_tree <- keep.tip(tree, lineageMbovis_tips)

# Save the subtrees as new nexus files
write.nexus(lineage1_tree, file='lineage1_tree.nexus')
write.nexus(lineage2_tree, file='lineage2_tree.nexus')
write.nexus(lineage3_tree, file='lineage3_tree.nexus')
write.nexus(lineage4_tree, file='lineage4_tree.nexus')
write.nexus(lineageMbovis_tree, file='lineageMbovis_tree.nexus')

ape::write.tree(lineage1_tree, file='lineage1_tree.newick')
ape::write.tree(lineage2_tree, file='lineage2_tree.newick')
ape::write.tree(lineage3_tree, file='lineage3_tree.newick')
ape::write.tree(lineage4_tree, file='lineage4_tree.newick')
ape::write.tree(lineageMbovis_tree, file='lineageMbovis_tree.newick')

# Make separate date files for each lineage as well

dates_all <- read.csv('Malawi_dates.csv')
 
dates_lineage1 <- dates_all[dates_all$name %in% lineage1_tips,]
dates_lineage2 <- dates_all[dates_all$name %in% lineage2_tips,]
dates_lineage3 <- dates_all[dates_all$name %in% lineage3_tips,]
dates_lineage4 <- dates_all[dates_all$name %in% lineage4_tips,]
dates_lineageMbovis <- dates_all[dates_all$name %in% lineageMbovis_tips,]

write.csv(dates_lineage1, 'dates_lineage1.csv')
write.csv(dates_lineage2, 'dates_lineage2.csv')
write.csv(dates_lineage3, 'dates_lineage3.csv')
write.csv(dates_lineage4, 'dates_lineage4.csv')
write.csv(dates_lineageMbovis, 'dates_lineageMbovis1.csv')











