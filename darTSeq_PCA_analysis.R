
# analysis of the SNPDart  using dartR 
#install.packages("devtools")

#install.packages("BiocManager")
#BiocManager::install(c("SNPRelate", "qvalue"))
#install_github("green-striped-gecko/dartR")
#install.packages("StAMPP")
#install.packages("Rtools")


#install.packages("PopGenome")
#install.packages("BiocManager")
#install.packages("vcfR")
#BiocManager::install("LEA")
#BiocManager::install("dartR")
#install.packages("pca3d")  # to plot 3D of the PCA


## load required libraries
library(devtools)
library(Rtools)
library(dartR)
library(StAMPP)  
library(pegas)
library(SNPRelate)
library(ggfortify)  # to autoplot PCA from Data
library(factoextra) # to plot PCA analysed
library(ape)
library(phangorn)
library(LEA)
library(pca3d)
library(MASS)

#install.packages("xfun")
#set workinfg dir

setwd("R:/CCDM-FGR-SE00542/Chala/hybrid/dart")


library(dartR)
# read the snp dataset
gl <- gl.read.dart(filename="snpDart.csv", topskip = 2,lastmetric ="Missing", ind.metafile = "snp.csv")##infact we can add su-population to see the variation explained by the differentiation of population 
gl



#generate reports
gl.report.reproducibility(gl)
# list the names
names(gl@other$loc.metrics)
head(gl@other$loc.metrics$Chr)
#listr individual names

indNames(gl)
write.csv(indNames(gl), file= "popName_group2.csv")

# plot the callRAte 

gl.report.callrate(gl)
## plot smear plot
gl.report.PICSnp(gl)
glPlot(gl)   
names(gl@other$loc.metrics)

# convert to plink format

gl2plink(gl, outfile = "plink2.txt", outpath = "R:/CCDM-FGR-SE00542/Chala/hybrid/dart", verbose = NULL)

#########################################FILE CONVERSION###################################################

gl2 <- gl.read.dart(filename="SNPDartFor_alder.csv", topskip = 2,lastmetric ="Missing", ind.metafile = "snp.csv") 
#convert to gds
gl2gds(gl2, outfile="dartR.gds",outpath = "R:/CCDM-FGR-SE00542/Chala/hybrid/dart")

dartgds <- snpgdsOpen("dartR.gds")

# convert to PED file of PLINK

snpgdsGDS2PED(dartgds,"dartSNPPED")

# convert to BED file of PLINK

snpgdsGDS2BED(dartgds,"dartSNP.bed")
# convert gds to EIGENSTRAT
snpgdsGDS2Eigen(dartgds,"dartSNP.eigen")
snpgdsClose(dartgds)
##############################################################
# write out the pop names to file

gl.make.recode.pop(gl, out.recode.file = "new_pop_assignments.csv",
                   outpath ="R:/CCDM-FGR-SE00542/Chala/hybrid/dart"
                   )
levels(pop(gl))
barplot(table(pop(gl)), las=2)
#############################################################################
# apply various filtering
# this is not applied in the current results
reproducible <- gl.filter.reproducibility(gl, t=99)
reproducible
monomorphic <- gl.filter.monomorphs(reproducible, v=0)
monomorphic
callRate <- gl.filter.callrate(monomorphic,method = "loc",threshold = 0.95)
filtered.gl  <- gl.filter.callrate(callRate, method="loc", threshold = 0.95)
filtered.gl

#reports after filtering 
gl.report.reproducibility(filtered.gl)
gl.report.callrate(filtered.gl)
gl.report.diversity(filtered.gl)
gl.report.monomorphs(filtered.gl)
gl.report.RepAvg(filtered.gl)

########################################################################
# check with wes as he said 5% missing were removed
callRate <- gl.filter.callrate(gl,method = "loc",threshold = 0.95)
dim(callRate)
#######################################################################

geneticDistance2 <- gl.dist.pop(filtered.gl, method="jaccard")
geneticDistance2

write.csv(as.data.frame(as.matrix(geneticDistance2)), file= "filtered_jaccard_genetic_distance.csv")

pwfst2 <-stamppFst(filtered.gl, nboots=100, percent=95, nclusters=8)  # nclusters is no of CPU to use

pwfst2
write.csv(pwfst2$Fsts, file = "filtered_genetic_distance_Fst_values.csv")

pwGst2 <-stamppNeisD(filtered.gl,pop=TRUE) ## pop=TRUE to compute at population level, FALSE at individual level
pwGst2

round(pwGst2,3)

# Calculate genetic distance between individuals 

Gent.dist2 <- stamppNeisD(filtered.gl, FALSE)  ## this copute for individuals as needed for AMOVA

Gent.dist2
round(Gent.dist2,3)

write.csv(as.data.frame(as.matrix(Gent.dist2)), file= "filtered_Nei_stance.csv")


#############################################################################

## genetic distance analysis
geneticDistance <- gl.dist.pop(gl, method="jaccard")
geneticDistance


write.csv(as.data.frame(as.matrix(geneticDistance)), file= "jaccard_genetic_distance.csv")

##  calculates pairwise Fst values between populations
# method 1 use stamppFst

pwfst <-stamppFst(gl, nboots=100, percent=95, nclusters=8)  # nclusters is no of CPU to use
pwfst
## method 2 use stamppNeisD (compute at individual or population level)
pwGst <-stamppNeisD(gl,pop=TRUE) ## pop=TRUE to compute at population level, FALSE at individual level
pwGst

round(pwGst,3)
##
pwGst <-stamppNeisD(gl,pop=TRUE)

# Calculate genetic distance between individuals 

Gent.dist <- stamppNeisD(gl, FALSE)  ## this copute for individuals as needed for AMOVA

Gent.dist
round(Gent.dist,3)

write.csv(as.data.frame(as.matrix(Gent.dist)), file= "Nei_stance.csv")


# Calculate AMOVA

amova_result <- stamppAmova(Gent.dist, gl, 100)
amova_result

head(amova_result$call)
head(amova_result$tab)
head(amova_result$varcoef)
head(amova_result$varcomp)

# write the result to file

write.csv(pwfst$Fsts, file = "genetic_distance_Fst_values.csv")
write.csv(pwfst$Pvalues, file = "genetic_distance_pvalues_values.csv")
write.csv(pwfst$Bootstraps, file = "genetic_distance_bootstrap100_values.csv")

# run PCA analysis

pc <- gl.pcoa(gl,nfactors=4)
pc


pc$eig
pc$call

names(pc)

#write the pc to fiel 
write.csv(pc$scores, file = "pca_scores.csv")
write.csv(pc$eig, file = "pca_eig.csv")
write.csv(pc$loadings, file = "pca_loding.csv")
write.csv(pc$call, file = "pca_CALL.csv")
# plot the fraction of variation explianed by  each PC
par(mfrow=c(2, 2))
barplot(pc$eig/sum(pc$eig)*100, ) 

## Plot populations PC # PCA individual as intities
gl.pcoa.plot(pc, gl, xaxis=1, yaxis=2) # labels="pop"
par(mfrow=c(2, 2))
pcoa
#plot color coded population
gl.pcoa.plot(pc, gl, labels="interactive", xaxis=1, yaxis=2)
ggplotly()

#screen plot
gl.pcoa.scree(pcoa) 

## create 3D  view of the file 

gl.pcoa.plot.3d(pc, gl)


library(pca3d)

# check test data 

data(metabo)
pca <- prcomp(metabo[,-1], scale.=TRUE)
gr <- factor(metabo[,1])
pca3d(pca, group=gr)
snapshotPCA3d(file="first_plot.png")
head(pca)

############MY data 

my4pca <- pc$scores  # copy the PCOAs 
my4pca
par(mfrow=c(2, 2))

group <- read.csv("snp.csv")
head(group)

pca3d(my4pca, radius = 5,show.shapes=TRUE,
      group = group$pop,legend="right",
      show.group.labels=FALSE,
      show.plane=FALSE,
      shape = 2,
      col= group$colour
      )


snapshotPCA3d(file="snp3D_corrected.png")


## compare with whole genome based snp analysis
# read vcf file from maf 

vcf.fn <- "R:/CCDM-FGR-SE00542/Chala/hybrid/dart/VCF/all_chromLevelSNP_rename_to_Chromosome.vcf"
head(vcf.fn)
# Reformat
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only")

snpgdsSummary("test3.gds")

set.seed(1000)
genofile <- snpgdsOpen("test3.gds")
pca <- snpgdsPCA(genofile,num.thread=8,autosome.only=FALSE)
pca
plot(pca$eigenvect[,1],pca$eigenvect[,2] ,col=as.numeric(substr(pca$sample, 1,3) == 'CCM')+3, pch=2)



#convert genom vcf to  bed file 
# this will write to file 
#snpgdsGDS2PED(genofile, ped.fn="genomeVcfped.ped")


# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
tab
# Draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")


# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.id
pop_code <- scan("R:/CCDM-FGR-SE00542/Chala/hybrid/dart/VCF/pop.txt", what=character())
pop_code
head(cbind(sample.id, pop_code))
write.csv((cbind(sample.id, pop_code)), file= "popName_group.csv")
# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
# Draw #############################################################
plot( tab$EV1,tab$EV2, col=as.integer(tab$pop), 
      xlab="PCoA axis 1 (79%)",cex = 3,
     ylab="PCoA Axis 2 (4%)", pch=17,
     axes = TRUE,
     )
abline(h=0, col="black")
abline(v=0, col="black")
legend("bottomright", legend=levels(tab$pop),  pch=19, col=1:nlevels(tab$pop))
####################################################################


# or in 3D plot 
# create tabe
tab3d <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  PCA1 = pca$eigenvect[,1],    # the first eigenvector
                  PCA2 = pca$eigenvect[,2],
                  PCA3 = pca$eigenvect[,3],  # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab3d)
write.csv(tab3d, file = "genomePCA.csv",)

# convert df to matrix

xx <- data.matrix(tab3d)
# get pop info

mygeopca <- read.csv("genomePCA.csv")
head(mygeopca)
mygeopca
colours  <- c("black","black", "red","black", "green","green","green","green","green","green","green","green")
pca3d(xx, radius = 4,show.shapes=TRUE,
      show.group.labels=FALSE,
      group =mygeopca$pop, 
      show.plane=FALSE,
      shape = 2,
      legend="right",
      col=colours,
      )
  
snapshotPCA3d(file="genomeScale_edited_3D.png")

xx

write.csv(tab3d, file = "genomePCA.csv",)


#play with parameters
pca3d(xx, radius = 3,show.shapes=TRUE,
      show.group.labels=FALSE,
      group =mygeopca$pop, 
      show.plane=FALSE,
      shape = 2,
      legend="right"
      )


snapshotPCA3d(file="genomeScale_3D.svg")






# define groups

par(mfrow=c(2, 2))
groups <- list(Ptt = sample.id[pop_code == "Ptt"],
               ptm = sample.id[pop_code == "Ptm"],
               Ptm_HR = sample.id[pop_code == "Ptm_HR"])

groups

prop <- snpgdsAdmixProp(pca, groups=groups)

# draw
plot(prop[, "Ptt"], prop[, "Ptm"], col=as.integer(tab$pop),
     xlab = "Admixture Proportion from Ptt",
     ylab = "Admixture Proportion from Ptm")
abline(v=0, col="gray25", lty=2)
abline(h=0, col="gray25", lty=2)
abline(a=0.5, b=-1, col="gray25", lty=2)
legend("topright", legend=levels(tab$pop), pch="o", col=1:4)
# draw
snpgdsAdmixPlot(prop, group=pop_code)





RV <- snpgdsEIGMIX(genofile)

snpgdsClose(genofile)

#example

# run eigen-analysis
RV <- snpgdsEIGMIX(genofile)








gl.tree.nj(gl, type="fan")  ## neibour joining tree

#CONVERT to GENOTYPE DATA

genofile <- snpgdsOpen("test3.gds")

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))
# run eigen-analysis
eig <- snpgdsEIGMIX(genofile)

# close the file
snpgdsClose(genofile)

#mygendata <- snpgdsGetGeno(genofile )
#head(mygendata)









snpgdsClose(genofile)
## analysis completed!







#filter the data 


gl2 <- gl.filter.repavg(gl, t=0.5)





#convert to newHtbrid testing file Format

glnew <- gl.nhybrids(gl)


# glnew <- gl.edit.recode.pop(gl)  # to edit the data 


###AMOVA EXAMPLES

# import genotype data and convert to allele frequecies 
data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, "r")

# Calculate genetic distance between individuals
potato.D.ind <- stamppNeisD(potato.freq, FALSE)

# Calculate AMOVA
stamppAmova(potato.D.ind, potato.freq, 100)


potato.freq


#chro1

myfas <- read.FASTA("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/Chr01_concatenated01 - Copy.fasta", type="DNA")

chr01.dist <- dist.dna(myfas, model='raw')

write.csv(as.data.frame(as.matrix(chr01.dist)), file= "chr01_distance.csv")

plot(chr01.dist, type = "phylogram", use.edge.length = TRUE,
      node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
      edge.color = "black", edge.width = 1, edge.lty = 1, font = 3,
      cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE,
      root.edge = FALSE, label.offset = 0, underscore = FALSE,
      x.lim = NULL, y.lim = NULL, direction = "rightwards",
      lab4ut = NULL, tip.color = "black", plot = TRUE,
      rotate.tree = 0, open.angle = 0, node.depth = 1,
      align.tip.label = FALSE)


myfascho2 <- read.FASTA("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/Chr02_concatenated01 - Copy.fasta", type="DNA")
chr02.dist <- dist.dna(myfascho2, model='raw')
write.csv(as.data.frame(as.matrix(chr02.dist)), file= "chr02_distance.csv")
bionJTree <- bionj(chr02.dist) 
plot(bionJTree)  # plot the tree 

#to bootstrap
func <- function(x) bionj(x, "raw")

dat <- as.matrix(chr02.dist)
dat
bstree <-  bionj(chr02.dist) 

mybootstr <- boot.phylo(dat,chr02.dist, function(xx) bionj(bstree),B= 10, mc.core=8)

## boot.phylo(bionJTree, FUN = function(xx)  bionj(xx), trees =TRUE )

## how to boot straop tree ?? 

myfascho3 <- read.FASTA("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/Chr03_concatenated008 - Copy.fasta", type="DNA")


chr03.dist <- dist.dna(myfascho3, model='raw')

write.csv(as.data.frame(as.matrix(chr03.dist)), file= "chr03_distance.csv")

## tree example



cat("(((Strix_aluco:4.2,Asio_otus:4.2):3.1,",
    "Athene_noctua:7.3):6.3,Tyto_alba:13.5);",
    file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
plot(tree.owls)
unlink("ex.tre") # delete the file "ex.tre"


### Reead DNA alignment infasta format as matrix 

chr03 <- read.dna("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/Chr03_concatenated008 - Copy.fasta",format="fasta", as.matrix = TRUE)

#compute distance  and tree
mydist <- dist.dna(chr03, model ="raw")  # K80 is default
tre  < bionj(mydist)
tre
bb <- boot.phylo(tre, chr03, function(xx) nj(dist.dna(xx)), B=100,mc.cores = 1)
plot(tre)
nodelabels(bb)

dm <- dist.ml(chr03) 
treeNJ <- bionj(dm)
bb <- boot.phylo(treeNJ, chr03, function(xx) nj(dist.dna(xx)), B=100,mc.cores = 1)
plot(treeNJ)
nodelabels(bb)
# add bootstrap to phylo object

treeNJ$node.label <- bb
# write tree to file 
write.tree(treeNJ,file = "chr03_bs_100_binj_mldist.dnd", digits = 2)

#####################################END for bionj tree construction 

### compute ML tree, we need to reload the data,
# read the file as pydat class
mypydat <- read.phyDat("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/Chr03_concatenated008 - Copy.fasta",format="fasta")

dm <- dist.ml(mypydat)
treeNJ <- NJ(dm)
# fit ML model
fit = pml(treeNJ, data=mypydat)
fit
methods(class="pml")
# fit the model to default ( JC model by default)
fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)

bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE,
                   control = pml.control(trace = 0))
# plot the trees 

par(mfrow=c(2,1))
par(mar=c(1,1,3,1))
#plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
plotBS(fitJC$tree, bs, p = 70, type="p")
cnet <- consensusNet(bs, p=0.2)
plot(cnet, "2D", show.edge.label=TRUE)
title("b)")

# update to change model parameters to GTR 

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

##############################

##try line plot for Phi rate test

chro12_profile <- read.csv("R:/CCDM-FGR-SE00542/Chala/NIMBUS/sibeliaz/all_sibeliaz_out_k11/all_aln/Final/split/PhiTest/chr12_profile_subset.csv", header = FALSE)
par(mar=c(2,2,3,2))
plot(chro12_profile$V1,chro12_profile$V2, type = "l", col = "red", xlab = "position", ylab = "Rate", 
     main = "Phi rate on Chr012")



###############################
# ANALYSIS AT INDIVIDUAL LEVEL NEED TO REDEFINE THE POPULATION
#

gl_Individual  <- gl #copy and store the original dataset in glind

pop(gl_Individual) <- indNames(gl_Individual)# redefine the population information

individual_dist <- gl.dist.pop(gl_Individual, method="jaccard")

write.csv(as.data.frame(as.matrix(individual_dist)), file= "indivudual_jaccard_distance.csv")


## convert vcf to other format 
library(vcfR)
remove(mygenome_vcf)
mygenome.vcf <- read.vcfR("R:/CCDM-FGR-SE00542/Chala/hybrid/dart/VCF/all_chromLevelSNP_rename_to_Chromosome.vcf")

head(mygenome.vcf)
mygeno <- vcf2geno("mygenome.vcf")

## oops this does not work 
# see https://www.biostars.org/p/231892/

mygenome.ped <- read.csv("R:/CCDM-FGR-SE00542/Chala/hybrid/dart/VCF/all_chromLevelSNP_rename_to_Chromosome.ped", sep ='\t')

genoOut = ped2geno(mygenome.ped, output.file = mygenome.geno, force = TRUE)

