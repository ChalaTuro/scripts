BiocManager::install("Sushi")
#version
library("Sushi")
library(tidyverse)
library(dplyr)
library(ggplot2)

#devtools::install_github("dvera/conifur")
#devtools::install_github("dvera/gyro")
#devtools::install_github("dvera/converge")
#devtools::install_github("dvera/travis")
library(devtools)
library("conifur")
library("converge")
library("gyro")

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

#BiocManager::install("NanoMethViz")

library(NanoMethViz)



#install_version(NanoMethViz) 

#example
library('Sushi')
Sushi_data = data(package = 'Sushi')
data(list = Sushi_data$results[,3])

head(Sushi_DNaseI.bedgraph)
chrom = "chr11"
chromstart = 1650000
chromend = 2350000
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,colorbycol= SushiColors(5))

## ACTUAL DATA  CHECKING 
teb040 <- read.csv2("teb_example.bedgraph",
                   sep="\t", header = TRUE)

head(teb040)
dim(teb040)
# convert the values to numeric and integers

teb040$start <- as.integer(teb040$start)
teb040$end <- as.integer(teb040$end)
teb040$value <- as.numeric(teb040$value)
head(teb040)
# subset scaf_1
#scaf_1 <- subset(teb040,scaf_1$chrom =="chr1") 
#head(scaf_1)
#scaf_1 <- as.data.frame(scaf_1)
#head(scaf_1)

# DEFINE REGION OF INTEREST TO PLOT  ON CHR01
chrom = "chr1"
chromstart = 2800000
chromend = 3063935

plotBedgraph(teb040,
             chrom,
             chromstart,
             chromend,
             color= "black",
            #ymax = 1.05,
             range=c(0,1.0),
             )

labelgenome(chrom,chromstart,chromend,n=6,scale="Mb")
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=1)
axis(side=2,las=2,tcl=.5)
# above is my default

###### test ggplot 
head(teb040)
ggplot(teb040) + geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=value))



# subset region of interest  

chrom = "chr1"
chromstart = 2542417
chromend = 4500609
plotBedgraph(teb040,
             chrom,
             chromstart,
             chromend,
             colorbycol = SushiColors(5),
             ymax = 1.05
             )
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")
axis(side=2,las=2,tcl=.5)

# zoom into a region of interest

# arrange the plot 
layout(matrix(c(1,1,2,3),2, 2, byrow = TRUE))
par(mar=c(3,4,1,1))

# get zoom region
zoomregion1 = c(250961,293542)
zoomregion2 = c(606105,700426)


zoomsregion(zoomregion1,extend=c(0.01,0.13),wideextend=0.05,
            offsets=c(0,0.580))
zoomsregion(zoomregion2,extend=c(0.01,0.13),wideextend=0.05,
            offsets=c(0,0.580))

# plot zoomed region 

plotBedgraph(teb040,chrom,chromstart=zoomregion1[1],
               chromend=zoomregion1[2],colorbycol= SushiColors(5))
labelgenome(chrom,chromstart=zoomregion1[1],chromend=zoomregion1[2],
              n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
zoombox()
mtext("metylation",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

#########################################################
# Plot  the methylation for 040

# read water treated one
# convert the methylation frequency to bedgraph to speed up loading data 
water040_bg <- read.csv2("water_18FRG040.result.tsv.plus.freq.bedgraph2",
                         sep="\t", header = TRUE)
head(water040_bg)
# rename column
names(water040_bg)[names(water040_bg) == "chromosome"] <- "chrom"
names(water040_bg)[names(water040_bg) == "methylated_frequency"] <- "value"
water040_bg$start <- as.integer(water040_bg$start)
water040_bg$end <- as.integer(water040_bg$end)
water040_bg$value <- as.numeric(water040_bg$value)
head(water040_bg)
# filter ourt regions with less than 0.05
water040_bg_filtered <- subset(water040_bg, value > 0.00) 

hist(water040_bg$value)

#read in Ethoh trated 

etoh040_bg <- read.csv2("etoh_18FRG040.result.tsv.plus.freq.bedgraph",
                         sep="\t", header = TRUE)
head(etoh040_bg)


# rename column
names(etoh040_bg)[names(etoh040_bg) == "chromosome"] <- "chrom"
names(etoh040_bg)[names(etoh040_bg) == "methylated_frequency"] <- "value"
etoh040_bg$start <- as.integer(etoh040_bg$start)
etoh040_bg$end <- as.integer(etoh040_bg$end)
etoh040_bg$value <- as.numeric(etoh040_bg$value)
head(etoh040_bg)

hist(etoh040_bg$value)
etoh040_bg_filtered <- subset(etoh040_bg, value > 0.00) 

#read in teb treated 
teb040_bg <- read.csv2("teb_18FRG040.result.tsv.plus.bedgraph",
                      sep="\t", header = TRUE)
head(teb040_bg)
# renqme column
names(teb040_bg)[names(teb040_bg) == "chromosome"] <- "chrom"
names(teb040_bg)[names(teb040_bg) == "methylated_frequency"] <- "value"
teb040_bg$start <- as.integer(teb040_bg$start)
teb040_bg$end <- as.integer(teb040_bg$end)
teb040_bg$value <- as.numeric(teb040_bg$value)
head(teb040_bg)
teb040_bg_filtered <- subset(teb040_bg, value > 0.0) 
dim(teb040_bg_filtered)
# plot histogram of the teb metylation
hist(teb040_bg$value)

# read in gene
ptt040_gene <- read.csv2("18FRG040_18FRG040_augustus.hints.genes.Count",
                           sep="\t", header = TRUE)
#colnames(ptt040_gene) <- c('chrom', 'start', 'end', 'value')
head(ptt040_gene)
#ptt040_gene$start <- as.integer(ptt040_gene$start)
#ptt040_gene$end <- as.integer(ptt040_gene$end)
#ptt040_gene$value <- as.numeric(ptt040_gene$value)
#head(ptt040_gene)
#read in repeats
 
#get repeat data 
ptt040_repeat <- read.csv2("18FRG040_repeats_filtered.bed.count",
                           sep="\t", header = TRUE)
#colnames(ptt040_repeat) <- c('chrom', 'start', 'end', 'value')
head(ptt040_repeat)
#ptt040_repeat$start <- as.integer(ptt040_repeat$start)
#ptt040_repeat$end <- as.integer(ptt040_repeat$end)
#ptt040_repeat$value <- as.numeric(ptt040_repeat$value)
#head(ptt040_repeat)

plotBedgraph(ptt040_gene,
             chrom,
             chromstart,
             chromend,
             range=c(0,20),
             # ymax = 1.05,
             transparency=0.5,
             color='blue',
             #border='red'
             )

plotBedgraph(ptt040_repeat,
             chrom,
             chromstart,
             chromend,
             range=c(0,20),
             # ymax = 1.05,
             #transparency=0.750,
             color='red',
             transparency=0.5,
             #color=SushiColors(2)(2)[1],
             #border='blue',
             overlay=TRUE,
             rescaleoverlay=TRUE
             )



###########################PLOT ALL over lay  for CYP scaf ##########
# scaffold_8 is cyp51 region 

layout(matrix(c(1,1,1,1),2, 2, byrow = TRUE))
par(mar = c(3, 5, 2, 2))  # define margins of the graph

par(mfrow=c(2, 1))
chrom = "scaffold_8"  
chromstart = 1
chromend = 3193661

plotBedgraph(teb040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             # ymax = 1.05,
             transparency=0.1,
             color='Orange',
             #border='red'
             )
plotBedgraph(water040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             # ymax = 1.05,
             color='blue',
             transparency=0.50,
             #color=SushiColors(2)(2)[1],
            #border='blue',
             overlay=TRUE,
             rescaleoverlay=TRUE
)
plotBedgraph(etoh040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
            # ymax = 1.05,
             transparency=0.1,
             #color=SushiColors(2)(2)[2],
             color='purple',
             border='purple',
             overlay=TRUE,
             rescaleoverlay=TRUE
                )

# ADD REPEAT DENSITY 
#plotBedgraph(teb040_repeat,
#             chrom,
#             chromstart,
#             chromend,
#             range=c(0,1),
#             # ymax = 1.05,
#             #transparency=.95,
#             #color=SushiColors(2)(2)[2],
#             color='green',
#             border='black',
#             overlay=TRUE,
#             rescaleoverlay=TRUE
#)

#ADD CDS density 
#plotBedgraph(ptt040_CDS,
#             chrom,
##             chromstart,
#             chromend,
##             range=c(0,1),
#             # ymax = 1.05,
#             #transparency=.95,
#             #color=SushiColors(2)(2)[2],
#             color='black',
#             overlay=TRUE,
#             rescaleoverlay=TRUE
#            )
             
# add labes to the plot 
labelgenome(chrom,chromstart,
            chromend,n=5,scale="Mb",
            chromfont = 1)
mtext("Methylation frequency", side=2,line=1.75,cex=1,font=0.25)
axis(side=2,las=2,tcl=.5)
# add legend to the plot 
legend("top",inset=0.025,legend=c("Teb","H2O", "EtOH"),
       fill=c("black","blue","purple"),text.font=1,
       cex=1)






######## apparantly the opposite holds

par(mfrow=c(3, 1))

plotBedgraph(water040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             # ymax = 1.05,
             #transparency=.50,
             color='blue',
             #color=SushiColors(2)(2)[1],
)
labelgenome(chrom,chromstart,
            chromend,n=5,scale="Mb",
            chromfont = 1)
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=1)
axis(side=2,las=2,tcl=.5)

plotBedgraph(etoh040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             #ymax = 1.05,
             #transparency=.50,
             #color=SushiColors(2)(2)[2],
             color='purple',
             
)

labelgenome(chrom,chromstart,
            chromend,n=5,scale="Mb",
            chromfont = 1)
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=1)
axis(side=2,las=2,tcl=.5)

plotBedgraph(teb040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             #ymax = 1.05,
             #transparency=0.50,
             #color=SushiColors(2)(2)[1])
             color='red',
)
labelgenome(chrom,chromstart,
            chromend,n=5,scale="Mb",
            chromfont = 1)
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=1)
axis(side=2,las=2,tcl=.5)




## here the reason for EToh looks high but actually not high.





#########################################################

# filp and plot  DID  worked well at this stage 
par(mfrow=c(1, 1))
par(mfrow=c(2,1),mar=c(1,4,1,1))
# plot chr08
chrom = "scaffold_8"  
chromstart = 1
#chromend = 3193661
chromend = 1000000

plotBedgraph(etoh040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             range=c(0,1.05),
             color='purple')
axis(side=2,las=2,tcl=.2)
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=0.5)
legend("topleft",inset=0.025,legend=c("EtOH", "H20"),
       fill=c("purple","blue"),text.font=1,
       cex=1.0)
plotBedgraph(water040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             flip=TRUE,
             color='blue')
labelgenome(chrom,chromstart,
            chromend,side=3, n=5,scale="Mb",
            chromfont = 1)
axis(side=2,las=2,tcl=.25,at=pretty(par("yaxp")[c(1,2)]),
     labels=-1*pretty(par("yaxp")[c(1,2)]))

mtext("Metylation frequency",side=2,line=1.75,cex=1,font=0.5)


####################################FINAL PLOT FOR GENE repeat and Water 
pdf("040_plots.pdf")
plotBedgraph(ptt040_gene,
             chrom,
             chromstart,
             chromend,
             range=c(0,20),
             # ymax = 1.05,
             color = 'blue',
             transparency=0.5,
             #border='red'
)

axis(side=2,las=2,tcl=.2)
mtext("Density",side=2,line=1.75,cex=1,font=0.5)
legend("topright",inset=0.025,legend=c("Gene", "Repeat", "H2O"),
       fill=c("blue","red", 'black'),text.font=1,
       cex=1.0)
plotBedgraph(ptt040_repeat,
             chrom,
             chromstart,
             chromend,
             range=c(0,20),
             # ymax = 1.05,
             #transparency=0.750,
             color='red',
             transparency=0.25,
             #color=SushiColors(2)(2)[1],
             #border='blue',
             overlay=TRUE,
             rescaleoverlay=TRUE,
             )

plotBedgraph(water040_bg_filtered,
             chrom,
             chromstart,
             chromend,
             flip=TRUE,
             color='black',
             transparency=0.10,
             )

labelgenome(chrom,chromstart,
            chromend,side=3, n=5,scale="Mb",
            chromfont = 1)
axis(side=2,las=2,tcl=.25,at=pretty(par("yaxp")[c(1,2)]),
     labels=-1*pretty(par("yaxp")[c(1,2)]))
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=0.5)
dev.off()
##

head(water040_bg_filtered)
bgHist(water040_bg_filtered, chrom='scaffold_1',score= water040_bg_filtered$value)
water040_bg_filtered


table_ethoh <- file(package = "NanoMethViz", "subset_water_18FRG040.result.tsv")
head(table_ethoh)

plot_region(table_ethoh, "scaffold_8", 1, 3193661)


################ summurize freq  for 026

Teb_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/026/teb_17FRG026.result.tsv.freq",
                        sep="\t", header = TRUE,dec = ".", numerals = "no.loss")
dim(Teb_freq)
Teb_freq <- as.data.frame(Teb_freq)
Teb_freq$methylated_frequency  <- as.integer(Teb_freq$methylated_frequency) 

head(Teb_freq)
# below require deplyer 
hist(Teb_freq$called_sites_methylated)
sumary <- Teb_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                "methylated_frequency"), sum, na.rm = TRUE)
sumary


Etoh_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/026/Etoh_17FRG026.result.tsv.freq",
                      sep="\t", header = TRUE,dec = ".", numerals = "no.loss")
dim(Etoh_freq)

Etoh_freq <- as.data.frame(Etoh_freq)
Etoh_freq$methylated_frequency  <- as.integer(Etoh_freq$methylated_frequency) 
sumaryE <- Etoh_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                 "methylated_frequency"), sum, na.rm = TRUE)
sumaryE

Water_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/026/water_17FRG026.result.tsv.freq",
                       sep="\t", header = TRUE,dec = ".", numerals = "no.loss")

Water_freq <- as.data.frame(Water_freq)
dim(Water_freq)
Water_freq$methylated_frequency  <- as.integer(Water_freq$methylated_frequency) 
sumaryW <- Water_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                 "methylated_frequency"), sum, na.rm = TRUE)
sumaryW


df <- bind_rows(sumaryW,sumaryE, sumary)
df
df$treatment <- c('Water','Etoh','Teb')
df 
df <- as.data.frame(df)
df

#barplot(df$called_sites_methylated)


p<-ggplot(data=df, aes(x=treatment, y=called_sites_methylated ))+
  geom_bar(stat = "identity",width=0.25)
  
p

# keep the order of the data to be ploted instead of sorted by x-axis 
df$treatment <- factor(df$treatment, levels = df$treatment[order(df$called_sites_methylated)])
p1 <-ggplot(data=df, aes(x=treatment, y=called_sites_methylated ))+
  geom_bar(stat = "identity",width=0.25)

p1

#########plot fraction of methylated sites 

df$feq <- df$called_sites_methylated /df$called_sites*100
df

p<-ggplot(data=df, aes(x=treatment, y=feq ))+
  geom_bar(stat = "identity",width=0.25)

p




################ summurize freq  for 040


Teb04_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/040/teb_18FRG040.result.tsv.freq",
                      sep="\t", header = TRUE,dec = ".", numerals = "no.loss")
dim(Teb04_freq)
Teb04_freq <- as.data.frame(Teb04_freq)
Teb04_freq$methylated_frequency  <- as.numeric(Teb04_freq$methylated_frequency) 

head(Teb04_freq)
# below require deplyer 
#hist(Teb04_freq$called_sites_methylated)

sumary04 <- Teb04_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                 "methylated_frequency"), sum, na.rm = TRUE)
sumary04


# ETOH
Etoh040_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/040/etoh_18FRG040.result.tsv.freq",
                       sep="\t", header = TRUE,dec = ".", numerals = "no.loss")
dim(Etoh040_freq)

Etoh040_freq <- as.data.frame(Etoh040_freq)
Etoh040_freq$methylated_frequency  <- as.integer(Etoh040_freq$methylated_frequency) 
sumaryE04 <- Etoh040_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                 "methylated_frequency"), sum, na.rm = TRUE)
sumaryE04

Water040_freq <- read.csv2("R:CCDM-FGR-SE00542/Chala/METHYLOME_2021/040/water_18FRG040.result.tsv.freq",
                        sep="\t", header = TRUE,dec = ".", numerals = "no.loss")

Water040_freq <- as.data.frame(Water040_freq)
dim(Water040_freq)
Water040_freq$methylated_frequency  <- as.integer(Water040_freq$methylated_frequency) 
sumaryW04 <- Water040_freq %>%
  summarise_at(c("num_motifs_in_group", "called_sites", "called_sites_methylated",
                 "methylated_frequency"), sum, na.rm = TRUE)

sumaryW04

df2 <- bind_rows(sumaryW04,sumaryE04,sumary04 )
df2
df2$treatment <- c('Water','Etoh','Teb')
df2 
df2 <- as.data.frame(df2)
df2

#barplot(df$called_sites_methylated)

#df2 <- bind_rows(sumaryW04,sumaryE04, sumary04)
#df2
#df2$treatment <- c('Water','Etoh','Teb')
#df2 
#df2 <- as.data.frame(df2)
#df2
#barplot(df$called_sites_methylated)

p<-ggplot(data=df2, aes(x=treatment, y=called_sites_methylated ))+
  geom_bar(stat = "identity",width=0.25)

p

# keep the order of the data to be ploted instead of sorted by x-axis 
df2$treatment <- factor(df2$treatment, levels = df2$treatment[order(df2$called_sites_methylated)])

p2<-ggplot(data=df2, aes(x=treatment, y=called_sites_methylated ))+
  geom_bar(stat = "identity",width=0.25)

p2

df2$feq <- df2$called_sites_methylated /df2$called_sites*100
df2

p2<-ggplot(data=df2, aes(x=treatment, y=feq ))+
  geom_bar(stat = "identity",width=0.25)

p2



############plot freq IN NA window OF 1kb 

#read in teb treated 
teb040_MEAN <- read.csv2("teb_040Mean.bed",
                       sep="\t", header = FALSE)

head(teb040_MEAN)
names(teb040_MEAN)
# renqme column
names(teb040_MEAN)[names(teb040_MEAN) == "V1"] <- "chrom"
names(teb040_MEAN)[names(teb040_MEAN) == "V2"] <- "start"
names(teb040_MEAN)[names(teb040_MEAN) == "V3"] <- "end"
names(teb040_MEAN)[names(teb040_MEAN) == "V4"] <- "value"
head(teb040_MEAN)

teb040_MEAN$start <- as.integer(teb040_MEAN$start)
teb040_MEAN$end <- as.integer(teb040_MEAN$end)
teb040_MEAN$value <- as.numeric(teb040_MEAN$value)
head(teb040_MEAN)

##################PLOT THE MEAN freq in 1 KB
chrom = "scaffold_8"  
chromstart = 1
chromend = 3193661

plotBedgraph(teb040_MEAN,
             chrom,
             chromstart,
             chromend,
             range=c(0,1),
             color='red')
axis(side=2,las=2,tcl=.2)
mtext("Metylation frequency",side=2,line=1.75,cex=1,font=0.5)
legend("topleft",inset=0.025,legend=c("EtOH", "H20"),
       fill=c("purple","blue"),text.font=1,
       cex=1.0)

# DEFINE REGION OF INTEREST TO PLOT  ON CHR01


chrom = "scaffold_8"
chromstart = 2800000
chromend =   3193661

plotBedgraph(teb040_MEAN,
             chrom,
             chromstart,
             chromend,
             color= "black",
             ymax = 1.05,
             range=c(0,0.5),
)



