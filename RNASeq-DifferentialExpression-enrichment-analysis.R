
# note remember to start Rstudio as admin  
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("edgeR")  
#BiocManager::install("multtest")
BiocManager::install("topGO")
BiocManager::install("goseq")
BiocManager::install("pathview")

library(ggplot2)
library(edgeR)
library(topGO)
library(tidyverse)
library("GenomicFeatures")   
library(rtracklayer)
library(AnnotationHub) # not installed 
library(Rsamtools)  # ok 
library("goseq") # not installed 
library("GO.db")
library("knitr")
library("data.table")

library("network")
library("sna")
library("ggplot2")
library("rbokeh")
library("ggnetwork")
library("scales")
library("ggiraph")
library("DT")
library("clusterProfiler")
library("topGO")
library(KEGGgraph) # not installed  
library(readMappings) # not installed 
library("geneLenDataBase")  #ok
library("org.Dm.eg.db")  #not installed 
library("readxl")
library("pathview")  # not inastalled 
library(dplyr)
library(Rsubread)  # not installed 
library("xlsx")




## input data from readmapping 
#featureCounts -F GTF -g gene_id --ignoreDup -a ../Pyrenophora_teres_f.teres_K0103_gffread_converted.gtf \
#-o K0103FeatureCounts.txt -T 8  -G ../K0103_packbio_sorted.fasta  K0103CTL1.bam K0103CTL2.bam K0103CTL3.bam \
#K0103TRT1.bam K0103TRT2.bam K0103TRT3.bam


data1 = read.table('E:/DIFFEXPR/AlignmentBased/mergedBam/K0103FeatureCounts.txt', header=T, row.names=1, com='#',sep = "\t")

head(data1)
data = subset(data1, select = c(6,7,8,9,10,11) )
data
col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = factor(c(rep("Control", 3), rep("Treated", 3)))
conditions
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study
exp_study = calcNormFactors(exp_study)
exp_study = estimateDisp(exp_study)
exp_study
et=exactTest(exp_study, pair=c("Control", "Treated"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="Control", sampleB="Treated", result_table)
head(result_table)
#result_table$logFC = -1 * result_table$logFC
write.table(result_table, file='star2.Control_vs_Treated.edgeR.DE_results', sep='\t', quote=F, row.names=T)
write.table(rnaseqMatrix, file='star2.Control_vs_Treated.edgeR.count_matrix', sep='\t', quote=F, row.names=T)
# write to excel file format
write.xlsx(result_table,file='star3_Control_vs_Treated_edgeR_DE_results.xlsx',
           sheetName = "rnaMappingExpr", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

##  get the top 10 genes up regulated 

top20 <- topTags(et,n=20, adjust.method = "BH", sort.by = "PValue", p.value = 0.01 )
top20
write.xlsx(top20, file='K0103_topdiff_regulated_genes.xlsx',
           col.names = TRUE, row.names = TRUE, append = FALSE )



#####source("/home/ubuntu/app/trinityrnaseq-v2.11.0/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf("star.Control_vs_Treated.edgeR.DE_results.MA_n_Volcano2.pdf")
plot_MA_and_Volcano(rownames(result_table), result_table$logCPM, result_table$logFC, result_table$FDR)
dev.off()

plotMDS(exp_study)
plotMDS(data)
plotBCV(exp_study)
abline(h=c(-1, 1), col="blue")
et = exactTest(exp_study, pair=c("Control", "Treated"))
#result_table$logFC = -1 * result_table$logFC
head(result_table)
plotMD(et, main = " Treated vs Control in K0103")  # plot sig datas
# remember rownames are gene ids
#result_table = data.frame(sampleA="Control", sampleB="Treated", result_table)
result_table = data.frame(result_table)
head(result_table)
result_table
dim(result_table)
summary(decideTests(et))  # number of genes up/down at 5% FDR



###########################################end of diff expression

#subset only few genes for testing 
result_table2 <- subset(result_table, result_table$PValue <= 0.01 & result_table$FDR <= 0.01) 
dim(result_table2)
summary(decideTests(et,adjust.method = "BH",lfc=1, p.value = 0.05)) # BH is default
genselected <- rownames(result_table2)
genselected 
write.table(genselected, file='selectedgenes.txt', sep='\t', quote=F, row.names=F)
length(genselected)

DE <- read.table("E:/DIFFEXPR/AlignmentBased/k103/selectedgenes.txt", quote="\"", header=T)
str(DE$x)
DEgenes <-DE$x
DEgenes

#geneListA <- rownames(result_table) # needed for topGo construction
# below package needed to read gff or gtf file 
#BiocManager::install("GenomicFeatures")
#BiocManager::install("AnnotationHub")
# RNAseq GO annotation based on ref 
# get GO annotation for topGo analysis

mygeneID2GO <- readMappings('E:/DIFFEXPR/AlignmentBased/k103/K0103_topGO_Map.txt')

str(head(mygeneID2GO))
mygnames <- names(mygeneID2GO)
mydiffgens <- unlist(mygnames, use.names=FALSE)
mygnames
mygnames <- sub("-*$", "", mygnames)
mygnames
#geneListA <-  factor(as.integer(mygnames%in%result_table))
#str(geneListA)
## for topGo: all genes genes in test; geneSel = geneswith pvalue < 0.01

geneListA <- factor(as.integer(mygnames %in% DEgenes))
names(geneListA) <- mygnames
str(geneListA)
# Ontology names CC, MF, BP
GOdata <- new("topGOdata", ontology = "MF", 
              allGenes = geneListA,
              annot = annFUN.gene2GO, 
              gene2GO = mygeneID2GO)
GOdata

numGenes(GOdata)
a <- genes(GOdata) ## obtain the list of genes
a
sg <- sigGenes(GOdata) # list of significant genes
str(sg)
numSigGenes(GOdata)

graph(GOdata) ## returns the GO graph
ug <- usedGO(GOdata)
head(ug)
sel.terms <- sample(usedGO(GOdata), 10)
num.ann.genes <- countGenesInTerm(GOdata, sel.terms) ## the number of annotated genes
num.ann.genes
ann.genes <- genesInTerm(GOdata, sel.terms) ## get the annotations
head(ann.genes)
#The scores for all genes, possibly annotated
#with names of the genes, can be obtained using the method scoresInTerm().
ann.score <- scoresInTerm(GOdata, sel.terms)
head(ann.score)
ann.score <- scoresInTerm(GOdata, sel.terms, use.names = TRUE)
head(ann.score)

#some statistics for a set of GO terms
termStat(GOdata, sel.terms)
#Performing the test Fisher test
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

FisherResult2 <- GenTable(GOdata,Fisher = resultFisher,topNodes = 25)
FisherResult2

write.xlsx(FisherResult2,file='Top_25_nodesIndividualTerm_FisherTest2.xlsx',
           sheetName = "MF2", 
           col.names = TRUE, row.names = TRUE, append = T)

#e Kolmogorov-Smirnov (KS)

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)

KSResult2 <- GenTable(GOdata,KS = resultKS,topNodes = 25)
KSResult2

#elm test
test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.stat)

elmResult2 <-GenTable(GOdata,KS = resultElim,topNodes = 25)
elmResult2


## Summarising the results
allRes <- GenTable(GOdata, classic = resultFisher, KS = resultKS, elm = resultElim,
                   orderBy = "elm", ranksOf = "classic", topNodes = 25)

allRes

write.xlsx(allRes,file='Top3_25_nodes_test_Go_test_result.xlsx',
           sheetName = "BP", 
           col.names = TRUE, row.names = TRUE, append = T)

#Analysing individual GOs
goID <- allRes[10, "GO.ID"]
goID

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = "sampleFile", useInfo = "all", pdfSW = TRUE)



#############################################################################
#completed topGo analysis 
# run Goseq if necessary 
#####################################################################




############################GOSeq Analysis
#get de genes for Goseq
de.genes = as.integer(result_table$PValue[result_table$logFC!=0] < 0.05) 
#get anotaton file

annotationFile <- import("Pyrenophora_teres_f.teres_K0103.sorted.gtf","gtf")
annotationFile
txtdb <- makeTxDbFromGRanges(annotationFile)
transcipt_length <- transcriptLengths(txtdb, with.cds_len=T)
transcipt_length
gene_list <- genes(txtdb)
gene_list
gene_list <- as.data.frame(gene_list) # convert to df
dim(gene_list)

geneid <- rownames(result_table)
geneid#names(genes) = row.names(result_table$table[result_table$table$logFC != 0, ])  
de.genes

# match column of two tables and by gene id and get the seq length 

result_table$Leng <- gene_list$width[match(row.names(result_table), gene_list$gene_id)] 
Leng <- as.integer(result_table$Leng)

# get Go annotation for GOSeq analysis

godata <- read.table('D:/DIFFEXPR/AlignmentBased/k103/K0103_GOSeq.txt', sep="\t", header=T)
godata <- godata[, 1:2]
godata <-as.data.frame(godata)

dim(godata)

pwf = nullp(de.genes,bias.data=Leng)

#row.names(pwf) <- names(Leng)
#head(pwf)

#goseq(pwf,gene2cat = godata, method = "Wallenius" ) # IS  THE DEFAULT 

gosq_res = goseq(pwf,gene2cat = godata,method = "Wallenius")
gosq_res
enriched.GO = gosq_res$category[p.adjust(gosq_res$over_represented_pvalue,
                                         method = "fdr") < 5] 

head(enriched.GO)
#with no bias inclusion
GO_nobias = goseq(pwf, gene2cat = godata, method = "Wallenius")
head(GO_nobias)
head(GO_nobias)

enriched.GO2 = gosq_res$category[p.adjust(gosq_res$over_represented_pvalue,
                                           method = "BH") < 0.05] 
enriched.GO2

for(go	in	enriched.GO[1:5]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

# plot the GOseq
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2019/RNAseq/html/06_Gene_set_testing.html#goseq
gosq_res %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

#both results become almost same in alignment based analysis of GO term from interproscan

# end of the analysis 


