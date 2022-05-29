# download the packgae and install

git clone https://github.com/PhHermann/LDJump.git
R CMD build LDJump
R CMD INSTALL LDJump_<version>.tar.gz


#Example;

#https://github.com/PhHermann/LDJump/blob/master/Example/Example.R
#require(LDJump)
#results = LDJump("/pathToSample/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fa", alpha = 0.05, segLength = 1000,
#                 pathLDhat = "/pathToLDhat", pathPhi = "/pathToPhi", format = "fasta", refName = NULL)
#postscript("Results.eps", horiz = F)
#plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
#dev.off()
#
#
#
mv Ptt-Ptm_Cyp51A_Nucleotide_alignment_wes_16072019.fasta  CYP51_align.fasta 

# start analysis of recombination for all chromosomes in Ptt/Ptm 

## chr1

#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr01_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()


write.csv(results, file="ldjump_out.csv")


#2

#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr02_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()

write.csv(results, file="ldjump_out.csv")

#3
#!/usr/bin/env Rscript
require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr03_concatenated008.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()


write.csv(results, file="ldjump_out.csv")

#4
#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr04_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()

write.csv(results, file="ldjump_out.csv")

#5
#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr05_concatenated008.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()

write.csv(results, file="ldjump_out.csv")


#6
#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr06_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#7
#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr07_concatenated008.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#8



#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr08_concatenated008.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#9

#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr09_concatenated008.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#10

#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr10_concatenated008.fasta ",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#11

#!/usr/bin/env Rscript
require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr01_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")

#12

#!/usr/bin/env Rscript

require(LDJump)
results= LDJump("/home/ubuntu/APPS/ChrLevel/Chr12_concatenated01.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("Results.png")
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump")
dev.off()
write.csv(results, file="ldjump_out.csv")


#!/bin/bash

for file in chr*; 
	do mkdir $file.ld;
	cd $file.ld
	bash $file
done

## Analyze CYP51_align
 #!/usr/bin/env Rscript

require(LDJump)
CYp51 <- LDJump("/home/ubuntu/APPS/CYP51/CYP51_align.fasta",alpha = 0.05, segLength = 100,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
png("CYp51.png")

plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump: CYP51")
dev.off()
 
lapply(CYp51, function(x) write.table( data.frame(x), 'CYp51.csv'  , append= T, sep=',' ))





	
#!/usr/bin/env Rscript

require(LDJump)

#re-analysis  and plot improvements

results <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr01_concatenated01.fasta",alpha = 0.05, segLength = 5000,cores =12,pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results2  <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr02_concatenated01.fasta",
        alpha = 0.05, segLength = 5000,cores =12,
        pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results3 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr03_concatenated008.fasta",
        alpha = 0.05, segLength = 5000,cores =12,
        pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results4 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr04_concatenated01.fasta",
       alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",
       pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")


#plot result 

#png("chr1-4Final.png", width=1000, height=1100)
svg("chr1-4Final.svg")
par(mfrow=c(4,1))
#plot(results[[1]],col="red", xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump: Chr01")
plot(results[[1]],xlim = c(0,800), col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 01")
plot(results2[[1]],xlim = c(0,700),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 02")
plot(results3[[1]],xlim = c(0,800),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 03")
plot(results4[[1]],xlim = c(0,700),col="blue",frame.plot=F,xlab="Segments", ylab = "Estimated Recombination Rate", main = "Chromosome 04")

dev.off()

# export result
lapply(results, function(x) write.table( data.frame(x), 'chr1.csv'  , append= T, sep=',' ))
lapply(results2, function(x) write.table( data.frame(x), 'chr02.csv'  , append= T, sep=',' ))
lapply(results3, function(x) write.table( data.frame(x), 'chr03.csv'  , append= T, sep=',' ))
lapply(results4, function(x) write.table( data.frame(x), 'chr04.csv'  , append= T, sep=',' ))

########################


###########################
#analyze data 


require(LDJump)
results5 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr05_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
results6 <-  LDJump("/home/ubuntu/APPS/ChrLevel/Chr06_concatenated01.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
results7 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr07_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")
results8 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr08_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

#plot resulty

#png("chr05-08Final.png", width=800, height=1100)
svg("chr05-08Final.svg")
par(mfrow=c(4,1))
plot(results5[[1]],xlim = c(0,500),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 05")
plot(results6[[1]],xlim = c(0,500),col="blue",frame.plot=F,xlab=NA, ylab= "Estimated Recombination Rate", main = "Chromosome 06")
plot(results7[[1]],xlim = c(0,500),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 07")
plot(results8[[1]],xlim = c(0,300),col="blue",frame.plot=F,xlab="Segments", ylab = "Estimated Recombination Rate", main = "Chromosome 08")

dev.off()
#export result

lapply(results5, function(x) write.table( data.frame(x), 'chr5.csv'  ,append= T, sep=',' ))
lapply(results6, function(x) write.table( data.frame(x), 'chr06.csv'  ,append= T, sep=',' ))
lapply(results7, function(x) write.table( data.frame(x), 'chr07.csv'  ,append= T, sep=',' ))
lapply(results8, function(x) write.table( data.frame(x), 'chr08.csv'  ,append= T, sep=',' ))

## write the results to file
#Analyze 
#R
#
require(LDJump)
results9 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr09_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results10 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr10_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results11 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr11_concatenated008.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")

results12 <- LDJump("/home/ubuntu/APPS/ChrLevel/Chr12_concatenated01.fasta",alpha = 0.05, segLength = 5000,cores =12, pathLDhat="/home/ubuntu/APPS/LDhat",pathPhi ="/home/ubuntu/APPS/PhiPack/Phi",format = "fasta")


#plot results

svg("chr9_12Final.svg")
#png("chr9_12Final.png", width=1000, height=1100)
par(mfrow=c(4,1))
plot(results9[[1]], xlim = c(0,400),col= "blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main= "Chromosome 09")
plot(results10[[1]],xlim = c(0,400),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 10")
plot(results11[[1]],xlim = c(0,300),col="blue",frame.plot=F,xlab=NA, ylab = "Estimated Recombination Rate", main = "Chromosome 11")
plot(results12[[1]],xlim = c(0,200),col="blue",frame.plot=F,xlab="Segments", ylab = "Estimated Recombination Rate", main = "Chromosome 12")

dev.off()

#export results

lapply(results9, function(x) write.table( data.frame(x), 'chr5.csv'  , append= T, sep=',' ))
lapply(results10, function(x) write.table( data.frame(x), 'chr06.csv'  , append= T, sep=',' ))
lapply(results11, function(x) write.table( data.frame(x), 'chr07.csv'  , append= T, sep=',' ))
lapply(results12, function(x) write.table( data.frame(x), 'chr08.csv'  , append= T, sep=',' ))
svg("chr9_12Final.svg", width=1000, height=1100)
#png("chr9_12Final.png", width=1000, height=1100)
par(mfrow=c(4,1))

dev.off()

## analysis completed 









