#####################################
#            EVE MODEL              #
#####################################

#remotes::install_gitlab("sandve-lab/evemodel")http://127.0.0.1:26581/graphics/plot_zoom_png?width=1347&height=834
library(evemodel)
library(ape)
library(phytools)
library(dplyr)

#============#
# INPUT DATA #
#============#

# read tab separated table
exprTbl <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\Broad Institute\\Flexible Homeostasis\\RNA-seq\\Differential Expression\\9 species\\TPM_9-species_expinallspp.tsv")
foldchangeTbl <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\Broad Institute\\Flexible Homeostasis\\RNA-seq\\Differential Expression\\9 species\\Fold Change Score\\Gene-score_9-species_2ormoresp-41C.tsv")
foldchangeTbl[is.na(foldchangeTbl)] <- 0
head(exprTbl)
head(foldchangeTbl)

# The species phylogy is in Newick format which can be read using the read.tree function in the ape package:
speciesTree <- read.tree("C:\\Users\\lexlu\\OneDrive\\Biology\\Broad Institute\\Flexible Homeostasis\\RNA-seq\\Differential Expression\\9 species\\EVE\\9_species-list.nwk")

# plot the species tree
plot(speciesTree)
#add.scale.bar(x=1,y=2,length = 5)

################################
#                              #
#       Expression levels      #
#                              #
################################

# REMOVE WEIRD RAT AND SQUIRREL
exprTbl <- exprTbl[,-c(114,115,116,144,145,146)]

samples <- colnames(exprTbl)
samples <- gsub("_camel","-camel", samples)
samples <- gsub("_badger","-badger", samples)
species <- unlist(lapply(strsplit(samples,"[_]"),`[[`, 1))
individual <- unlist(lapply(strsplit(samples,"[_]"),`[[`, 2))
temperature <- unlist(lapply(strsplit(samples,"[_]"),`[[`, 3))
glucose <- unlist(lapply(strsplit(samples,"[_]"),`[[`, 4))

# REMOVE WEIRD RAT AND SQUIRREL
heat <- exprTbl[,which(temperature=="41C")]

samples_heat <- samples[which(temperature=="41C")]
species_heat <- species[which(temperature=="41C")]

# I can also remove genes for which species have 0 expression
# Calculate sum across each group for each row using rowsum
rowsum_heat <- rowsum(t(heat), group = species_heat)
has.expression_heat <- rowsum_heat>0
sum.has.expression_heat <- colSums(has.expression_heat)
# Select the genes with expression in all species
has.allspp.expression_heat <- which(sum.has.expression_heat==9)
all.spp.exp_heat <- heat[has.allspp.expression_heat,]

# The table needs to be converted to a matrix:

# The first column in the table is the ortholog group ID
# which will be the rownames of the matrix.
# TODO: The rownames are currently not returned in the results, but they should.
exprMat <- as.matrix(log2(heat+0.001))
exprMat.no.miss <- as.matrix(log2(all.spp.exp_heat+0.001))

#================================================================#
#                                                                #
#             Camels, 13lgs and bat are flexible                 #
#                                                                #
#================================================================#

#=================================================================#
# Mapping species in tree to the columns in the expression matrix #
#=================================================================#

# The evemodel methods needs to know which species each column in the expression matrix belongs to. This is done by creating a vector of species names corresponding to the tip labels in the species tree.

# the species names in the tree is given by the tip.labels
speciesTree$tip.label

# the columns in the expression matrix are:
colnames(exprMat)

# create a species vector that maps to the tree tip labels
species <- gsub("bactrian-camel","Camelus_bactrianus", species)
species <- gsub("bat","Myotis_lucifugus", species)
species <- gsub("dromedary-camel","Camelus_dromedarius", species)
species <- gsub("dolphin","Tursiops_truncatus", species)
species <- gsub("honey-badger","Mellivora_capensis", species)
species <- gsub("human","Homo_sapiens", species)
species <- gsub("rat","Rattus_norvegicus", species)
species <- gsub("rhino","Ceratotherium_simum", species)
species <- gsub("squirrel","Ictidomys_tridecemlineatus", species)
species <- species[which(temperature=="41C")]
species

#==============================#
# Running the theta shift test #
#==============================#

# Test dataset with only the first 100 genes
#exprMat.test <- exprMat[1:100,]

# We test two competing evolutionary hypotheses specified as two OU models, one with and one without a phenotype-associated theta shift.
# We then fit our data (gene expression from extant species) to the competing hypotheses of one shared theta (one theta for the entire gene tree, the gene tree describing the evolutionary relationship between genes in one ortholog group) and two thetas (??2 for the foreground clades and ??1 for the rest of the tree) using the evemodel

# select edges with theta shift
flexible <- c("Camelus_bactrianus","Camelus_dromedarius","Myotis_lucifugus","Ictidomys_tridecemlineatus")
flex_edges <- c(getEdgesFromMRCA(speciesTree,tips = c("Camelus_bactrianus","Camelus_dromedarius"),includeEdgeToMRCA = F),
                which.edge(speciesTree,"Myotis_lucifugus"),
                which.edge(speciesTree,"Ictidomys_tridecemlineatus"))
shiftyEdges = 1:Nedge(speciesTree) %in% flex_edges

plot(speciesTree,edge.color = ifelse(shiftyEdges,"red","black"))

# Running the two theta model
twoThetaRes <- twoThetaTest(tree = speciesTree, gene.data = exprMat,
                            isTheta2edge = shiftyEdges, colSpecies = species,
                            cores = 4)

# function to make the results in table format
twoThetaRes2table <- function(res, OGIDs){
  tibble( OG = OGIDs, 
          LRT = res$LRT) %>% 
    bind_cols(as_tibble(res$twoThetaRes$par)) %>% 
    rename( theta = theta1, thetaShift = theta2, sigma.sq = sigma2) %>% 
    mutate( shift.direction = ifelse(thetaShift>theta, "up","down")) %>% 
    mutate( pval = pchisq(LRT,df = 1,lower.tail = F)) %>% 
    mutate( padjust = p.adjust(pval,method = "fdr")) %>%
    mutate( llTwoTheta = res$twoThetaRes$ll ) %>% 
    mutate( llOneTheta = res$oneThetaRes$ll ) %>% 
    mutate( isSig = padjust < 0.05)
}

resTbl <- twoThetaRes2table(twoThetaRes,OGIDs = rownames(exprMat))
write.table(resTbl,"twoTheta.allflexlible-41C.tsv",sep = "\t",quote = F,row.names = F)

resTbl <- read.delim("twoTheta-41C.tsv",sep = "\t")

# Make plot
for(g in resTbl$OG[which(resTbl$pval<0.05)]){
  w <- which(rownames(exprMat)==g)
  xe<-setNames(exprMat[w,],species)
  pdf(paste0("significant genes\\",g,"_exp-41C.pdf"))
  par(mfrow=c(1,2))
  plotTree(speciesTree,mar=c(5.1,1.1,2.1,0.1))
  par(mar=c(5.1,0.1,2.1,1.1))
  boxplot(xe~factor(names(xe),levels=speciesTree$tip.label),horizontal=TRUE,
          axes=FALSE,xlim=c(1,Ntip(speciesTree)),col=c("lightblue","lightblue","red4","red4","lightblue","red4","red4","lightblue","lightblue"),
          boxwex=.5,xlab=NA)
  axis(1)
  title(xlab="log2 (gene expression at 41C)")
  title(g)
  dev.off()
}

# # attempt to plot icons for spps
# library(ggplot2)
# library(ggtree)
# library(ggimage)
# 
# sp.lt <- speciesTree$tip.label
# sp.lt[sp.lt=="Myotis_lucifugus"] <- "Lasiurus_cinereus"
# d <- ggimage::phylopic_uid(sp.lt)
# d$flexible <- c(30,30,80,80,30,80,80,30,30)
# 
# p <- ggtree(speciesTree) %<+% d + 
#   geom_tiplab(aes(image=uid, colour=flexible), geom="phylopic", offset=2.5) +
#   scale_color_viridis_c()


#==============================#
# Running the beta shared test #
#==============================#

sharedBetaRes <- betaSharedTest(tree = speciesTree, gene.data = exprMat.no.miss, colSpecies = species)

sharedBetaRes2table <- function(res, OGIDs){
  tibble( OG = OGIDs, 
          LRT = res$LRT) %>% 
    bind_cols(as_tibble(res$indivBetaRes$par)) %>%  
    mutate( shift.direction = ifelse(beta>res$sharedBeta, "relaxed","divergent")) %>% 
    mutate( pval = pchisq(LRT,df = 1,lower.tail = F)) %>%
    mutate( padjust = p.adjust(pval,method = "fdr")) %>%
    mutate( llindivBeta = res$indivBetaRes$ll ) %>% 
    mutate( llsharedBeta = res$sharedBetaRes$ll ) %>% 
    mutate( isSig = padjust < 0.05)
}

resBetaTbl <- sharedBetaRes2table(sharedBetaRes,OGIDs = rownames(all.spp.exp_heat))
sharedBetaRes$sharedBeta

write.table(resBetaTbl,"sharedBeta-allflexible-41C.tsv",sep = "\t",quote = F,row.names = F)

resBetaTbl <- read.delim("sharedBeta-41C.tsv",sep = "\t")

# Make plot
for(g in resBetaTbl$OG[which(resBetaTbl$pval<0.05)]){
  w <- which(rownames(exprMat.no.miss)==g)
  xe <- setNames(exprMat.no.miss[w,],species)
  pdf(paste0("significant genes\\",g,"beta-test_exp-41C.pdf"))
  par(mfrow=c(1,2))
  plotTree(speciesTree,mar=c(5.1,1.1,2.1,0.1))
  par(mar=c(5.1,0.1,2.1,1.1))
  boxplot(xe~factor(names(xe),levels=speciesTree$tip.label),horizontal=TRUE,
          axes=FALSE,xlim=c(1,Ntip(speciesTree)),col=c("lightblue","lightblue","red4","red4","lightblue","red4","red4","lightblue","lightblue"),
          boxwex=.5,xlab=NA)
  axis(1)
  title(xlab="log2 (gene expression at 41C)")
  title(g)
  dev.off()
}

###########################################
# Intersect two-theta and two-beta tests  #
###########################################

common <- intersect(resTbl$OG[resTbl$pval<0.05],resBetaTbl$OG[resBetaTbl$pval<0.05])
write.table(common,"Intersect_two-thetatwo-beta-41C_allflexible.tsv",sep = "\t")

#=================#
# Simulating data #
#=================#

# Given a set of parameters it is possible to simulate gene expression. But what are reasonable parameters?
# Let's see what the median parameters from the fitted results looks like:
apply(na.omit(twoThetaRes$oneThetaRes$par),2,median)

theta <- 3.95
sigma2 <- 0.24
alpha <- 0.04
betaShared <- 0.44

#set.seed(123) # we want the same simulation each time...

null.LRT <- numeric()
for(i in 1:100){
  print(paste(i,"of",100))
  nullData_1000genes <- simOneTheta(n = 1000, tree = speciesTree, colSpecies = species,
                                    theta = theta, sigma2 = sigma2, alpha = alpha, beta = betaShared)
  
  # The simulated data can then be used to run the tests again:
  test.nullData_full <- twoThetaTest(tree = speciesTree, gene.data = nullData_1000genes, 
                                     isTheta2edge = shiftyEdges, colSpecies = species, cores=4) 
  
  # Get LRT scores
  test.nullData_full.LRT <- test.nullData_full$LRT
  
  # Append to LRT vector
  null.LRT <- c(null.LRT,test.nullData_full.LRT)
  
  #test.nullData_betaShared <- betaSharedTest(tree = speciesTree, gene.data = nullData_1000genes,
  #                                           colSpecies = species, cores=4)
}

LRT.null.threshold <- quantile(na.omit(null.LRT),0.95)
p.null.threshold <- dchisq(LRT.null.threshold,df = 1)

################################
#                              #
#          Fold Change         #
#                              #
################################

# The table needs to be converted to a matrix:
fold.changeMat <- as.matrix(foldchangeTbl)

#================================================================#
#                                                                #
#             What is going on with flexible spp?                #
#                                                                #
#================================================================#

#=================================================================#
# Mapping species in tree to the columns in the expression matrix #
#=================================================================#

# The evemodel methods needs to know which species each column in the expression matrix belongs to. This is done by creating a vector of species names corresponding to the tip labels in the species tree.

# the species names in the tree is given by the tip.labels
speciesTree$tip.label

# the columns in the expression matrix are:
colnames(fold.changeMat)

# create a species vector that maps to the tree tip labels
species <- gsub("bactrian.camel","Camelus_bactrianus", species)
species <- gsub("bat","Myotis_lucifugus", species)
species <- gsub("dromedary.camel","Camelus_dromedarius", species)
species <- gsub("dolphin","Tursiops_truncatus", species)
species <- gsub("honey.badger","Mellivora_capensis", species)
species <- gsub("human","Homo_sapiens", species)
species <- gsub("rat","Rattus_norvegicus", species)
species <- gsub("rhino","Ceratotherium_simum", species)
species <- gsub("squirrel","Ictidomys_tridecemlineatus", species)
species <- species[which(temperature=="41C")]
species

#==============================#
# Running the theta shift test #
#==============================#

# We test two competing evolutionary hypotheses specified as two OU models, one with and one without a phenotype-associated theta shift.
# We then fit our data (gene expression from extant species) to the competing hypotheses of one shared theta (one theta for the entire gene tree, the gene tree describing the evolutionary relationship between genes in one ortholog group) and two thetas (??2 for the foreground clades and ??1 for the rest of the tree) using the evemodel

# select edges with theta shift
flexible <- c("Camelus_bactrianus","Camelus_dromedarius","Myotis_lucifugus","Ictidomys_tridecemlineatus")
flex_edges <- c(getEdgesFromMRCA(speciesTree,tips = c("Camelus_bactrianus","Camelus_dromedarius"),includeEdgeToMRCA = F),
                which.edge(speciesTree,"Myotis_lucifugus"),
                which.edge(speciesTree,"Ictidomys_tridecemlineatus"))
shiftyEdges = 1:Nedge(speciesTree) %in% flex_edges

# flexible <- c("Camelus_bactrianus","Camelus_dromedarius")
# flex_edges <- c(getEdgesFromMRCA(speciesTree,tips = c("Camelus_bactrianus","Camelus_dromedarius"),includeEdgeToMRCA = F))
# shiftyEdges = 1:Nedge(speciesTree) %in% flex_edges

plot(speciesTree,edge.color = ifelse(shiftyEdges,"red","black"))

# Running the two theta model
twoThetaRes_foldchange <- twoThetaTest(tree = speciesTree, gene.data = fold.changeMat,
                            isTheta2edge = shiftyEdges, colSpecies = species,
                            cores = 4)

# function to make the results in table format
twoThetaRes2table <- function(res, OGIDs){
  tibble( OG = OGIDs, 
          LRT = res$LRT) %>% 
    bind_cols(as_tibble(res$twoThetaRes$par)) %>% 
    rename( theta = theta1, thetaShift = theta2, sigma.sq = sigma2) %>% 
    mutate( shift.direction = ifelse(thetaShift>theta, "up","down")) %>% 
    mutate( pval = pchisq(LRT,df = 1,lower.tail = F)) %>% 
    mutate( padjust = p.adjust(pval,method = "fdr")) %>%
    mutate( llTwoTheta = res$twoThetaRes$ll ) %>% 
    mutate( llOneTheta = res$oneThetaRes$ll ) %>% 
    mutate( isSig = padjust < 0.05)
}

resTbl.fold_change <- twoThetaRes2table(twoThetaRes_foldchange,OGIDs = rownames(fold.changeMat))
write.table(resTbl.fold_change,"twoTheta-41C_foldchange.tsv",sep = "\t",quote = F,row.names = F)

# plot(fold.changeMat[which(rownames(fold.changeMat)=="ENSG00000113068"),],col=c(rep("black",5),
#                           rep("red",5),
#                           rep("blue",5),
#                           rep("yellow",5),
#                           rep("green",5),
#                           rep("pink",5),
#                           rep("orange",5),
#                           rep("purple",4),
#                           rep("brown",4)),
#      pch=c(rep(17,5),rep(16,5),rep(17,5),rep(16,28)))

for(g in resTbl.fold_change$OG[which(resTbl.fold_change$pval<0.05)]){
  w <- which(rownames(fold.changeMat)==g)
  xe<-setNames(fold.changeMat[w,],species)
  pdf(paste0("significant genes\\",g,"_foldchange-41C.pdf"))
  par(mfrow=c(1,2))
  plotTree(speciesTree,mar=c(5.1,1.1,2.1,0.1))
  par(mar=c(5.1,0.1,2.1,1.1))
  boxplot(xe~factor(names(xe),levels=speciesTree$tip.label),horizontal=TRUE,
          axes=FALSE,xlim=c(1,Ntip(speciesTree)),col=c("lightblue","lightblue","red4","red4","lightblue","red4","red4","lightblue","lightblue"),
          boxwex=.5,xlab=NA)
  axis(1)
  title(xlab="log2 (fold change at 41C)")
  title(g)
  dev.off()
}

# Let's make a heatmap of only the significant genes

sig.genes <- na.exclude(resTbl.fold_change$OG[resTbl.fold_change$pval<0.01])
sig.fold.changeMat <- fold.changeMat[row.names(fold.changeMat) %in% sig.genes,]
rowsum_heat <- aggregate(t(sig.fold.changeMat), by = list(species), FUN = median)

fold.change.41.sign <- data.frame(t(rowsum_heat[,-1]))
colnames(fold.change.41.sign) <- rowsum_heat$Group.1

fold.change.41.dist <- dist(as.matrix(t(fold.change.41.sign)))
plot(nj(fold.change.41.dist))

###################
# Produce heatmap #
###################

library(cluster)
library(gplots)

heatmap.2(as.matrix(t(fold.change.41.sign)),
          col = viridis::viridis_pal(),
          trace="none",
          scale="column",
          density.info="none",
          cexRow=0.75,
          cexCol=0.3,
          sepcolor = "white",
          dendrogram = "row",
          hclustfun = function(x) hclust(x, method="ward.D"),
          #colsep=1:ncol(t(fold.change.41.sign)), # Add vertical grid lines
          rowsep=1:nrow(t(fold.change.41.sign)))
