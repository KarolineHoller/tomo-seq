# Description ####
# How we analyzed tomo-seq data for zebrafish one-cell stage embryos (sample here is called tomo13)

#Versions:
#R version 3.2.3 (2015-12-10)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: OS X 10.14.6 (unknown)
#
#locale:
#  [1] C

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] gridExtra_2.3   kohonen_3.0.4   zoo_1.8-1       tsne_0.1-3      cluster_2.0.6   pheatmap_1.0.12
#[7] reshape2_1.4.3  ggplot2_2.2.1  

#loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.14       lattice_0.20-35    digest_0.6.13      MASS_7.3-50        grid_3.2.3        
#[6] plyr_1.8.4         gtable_0.2.0       magrittr_1.5       scales_0.5.0       rlang_0.1.6       
#[11] stringi_1.2.4      lazyeval_0.2.1     labeling_0.3       RColorBrewer_1.1-2 tools_3.2.3       
#[16] stringr_1.3.1      munsell_0.5.0      colorspace_1.3-2   tibble_1.3.4      

setwd("/your/file/path/tomoseq.zebrafish.rep1.csv")
list.files()
#The needed files for that script are count tables from the mapping script. I put these count tables in the upon #mentioned folder.

# Parameters ####
min.transcripts = 8000
min.e = 5
min.s = 1

# Dependencies ####
library(ggplot2)
library(reshape2)
library(pheatmap)
library(cluster)
library(tsne)
library(zoo)
library(kohonen)
library(gridExtra)
require(VennDiagram)


# Functions ####

plot.gene.name <- function(name, expression.data, exact = F){
  if(exact == T){
    gene.expression <- expression.data[expression.data$Gene == name, ]
    plot.title <- paste("Spatial expression of", gene.expression[1])
    plot.data <- 
      data.frame(Section = as.integer(sapply(colnames(gene.expression)[-1], function(x) unlist(strsplit(x, "X"))[2])),
                 Expression = t(gene.expression[-1]))
    colnames(plot.data)[2] <- "Expression"
    print(
      ggplot(plot.data, aes(x = Section, y = Expression)) +geom_bar(stat = "identity") +
        theme_light() +
        theme(text = element_text(size=15)) +
        labs(title = plot.title)
    )
  }else{
    gene.lines <- grep(name, expression.data$Gene)
    for(line in gene.lines){
      gene.expression <- expression.data[line, ]
      plot.title <- paste("Spatial expression of", gene.expression[1])
      plot.data <- 
        data.frame(Section = as.integer(sapply(colnames(gene.expression)[-1], function(x) unlist(strsplit(x, "X"))[2])),
                   Expression = t(gene.expression[-1]))
      colnames(plot.data)[2] <- "Expression"
      print(
        ggplot(plot.data, aes(x = Section, y = Expression)) + geom_bar(stat = "identity") +
          theme_light() + 
          theme(text = element_text(size=15)) +
          labs(title = plot.title)
      )
    }
  }
}

euc.dist <- function(gene, expression.data){
  i=grep(gene, expression.data$Gene)
  m=as.matrix(dist(rbind(expression.data[i,], expression.data)[,-1], upper=T))
  d=data.frame("Gene"=expression.data[rownames(m), "Gene"],
               "distance"=m[,1])
  d.s=d[order(d$distance, decreasing = F),]
  d.s=d.s[-2,]
  return(d.s)
}

rpm=function(lib.size, read.counts){
  read.n=read.counts/lib.size*10^6
  return(read.n)
}

smooth.expression=function(gene.expression.data, bin.length){
  f=c(sum(gene.expression.data[1:2])/2,"",sum(gene.expression.data[(length(gene.expression.data)-1):length(gene.expression.data)])/2)
  genes.s=t(rollmean(t(gene.expression.data), k=bin.length, align = 'center', fill=f))
  return(genes.s)
}

plot.gene <- function(name, expression.data){
    gene.expression <- expression.data[expression.data$Gene == name, ]
    plot.title <- paste("", gene.expression[1])
    plot.data <- 
      data.frame(Section = as.integer(sapply(colnames(gene.expression)[-1], function(x) unlist(strsplit(x, "X"))[2])),
                 Expression = t(gene.expression[-1]))
    colnames(plot.data)[2] <- "relative Expression"
    p=ggplot(plot.data, aes(x = Section, y = Expression)) +geom_bar(stat = "identity") +
        theme_light() +
        theme(text = element_text(size=15)) +
        labs(title = plot.title)
    return(p)
}

# Load data ####
transcripts.in.tomo13 <- read.table("/your/file/path/tomoseq.zebrafish.rep1.transcripts.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(transcripts.in.tomo13)[-1]=paste("X", 1:96, sep = "")
transcripts.tomo13 <- transcripts.in.tomo13[!grepl("ERCC", transcripts.in.tomo13$GENEID), ]
#counts.in.tomo13 <- read.table("/your/file/path/tomoseq.zebrafish.rep1.counts.csv", sep = "\t", header = T, stringsAsFactors = F)
#counts.tomo13 <- counts.in.tomo13[!grepl("ERCC", counts.in.tomo13$GENEID), ]
UMIs.in.tomo13 <- read.table("/your/file/path/tomoseq.zebrafish.rep1.UMIs.csv", sep = "\t", header = T, stringsAsFactors = F)
UMIs.tomo13 <- UMIs.in.tomo13[!grepl("ERCC", UMIs.in.tomo13$GENEID), ]

add.names= read.table("/your/file/path/gene_names_translator.csv", sep = ",", header=T, stringsAsFactors = F) #that file translates GENEIDs to gene names, downloaded from ensemble, matching the transcriptome that we mapped to

#table 1 contains normalized expression and SOM profile per gene for three different samples of one-cell stage zebrafish embryos:
tomo.replicates.SOM <- read.table("/your/file/path/table1.csv", sep = ",", header = T, stringsAsFactors = F)

# Basic inspection of data quality ####
#ERCC content per section in tomo13 - checks if the internal controls are distributed evenly (only few dropouts are a mark of high quality)
ERCC.tomo13.s=colSums(ERCC.tomo13[,-1])
#this produces the plot in Figure S1c:
plot(ERCC.tomo13.s, type="h",
     xlab="Section number", ylab="ERCC reads",
     main="ERCC reads in tomo13")

# Sequencing statistics tomo13####
#the number of reads per UMI are a measure of how deeply a library is sequenced - desirable would be that every UMI has been seen multiple times
#reproduces plot in Fig. 1d
df_overseq <- data.frame("gene" = counts.tomo13$GENEID,
                         "reads" = rowSums(counts.tomo13[,-1]),
                         "UMI" = rowSums(UMIs.tomo13[,-1]))

#function for UMI saturation (Gruen et al., 2014):
#while 4^6*96 is the maximum amount of UMIS that we could potentially recover, with 6 bp random
#barcode, multiplied with 96 sections
UMI.saturation <- function(x) {
  (1- exp(-x/(4^6*96)))*4^6*96
}

#this code produces the plot in Figure 1d
ggplot(data = df_overseq) +
  theme_light() +
  scale_x_log10(limits = c(1, 10^8)) +
  scale_y_log10() +
  geom_point(aes(x= reads, y=UMI), alpha = 0.3, size= 0.5) +
stat_function(fun = UMI.saturation, col = "turquoise4") + #maximum complexity
geom_hline(yintercept= 1, col="magenta4")  #minimum complexity
#ggsave("tomo13.UMI.saturation.model.png", width = 3.5, height = 2.5, unit="in")

bc.stats.tomo13 <- data.frame(Section = factor(c(1:96)),
                       Transcripts = colSums(transcripts.tomo13[, -1]),
                       Genes = apply(transcripts.tomo13[, -1], 2, function(x) sum(x > 0)),
                       ERCC=ERCC.tomo13.s)

#visualizes the distribution of transcript- and read recovery per section as shown in Figure 1b
ggplot(bc.stats.tomo13) + 
  theme_bw() +
  geom_line(aes(x = Section, y = Transcripts, group=1)) +
  #scale_y_log10() + 
  geom_hline(yintercept=min.transcripts, color="red") +
  geom_line(aes(x=Section, y=ERCC, group=2), color="blue")+
  theme(axis.text.x = element_text(angle = 45)) +
  labs(title = "Transcript recovery in tomo13")
#ggsave("transcripts_per_section_tomo13.logcut8000.pdf", width = 12, height = 6, unit="in")

ggplot(bc.stats.tomo13, aes(x = Section, y = Genes, group=1)) + geom_line() +
  #scale_y_log10() + 
  #geom_hline(yintercept=min.transcripts, colour="red") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(title = "Gene recovery in tomo13")
#ggsave("genes_per_section_tomo13.png", width = 12, height = 6, unit="in")

#Plot the reads per section after filtering out low-qualitly sections
bc.stats.2.tomo13 <- bc.stats.tomo13
bc.stats.2.tomo13 <- bc.stats.tomo13[bc.stats.tomo13$Transcripts > min.transcripts, ]
#that filters out sections with low quality, e.g. low recovery of transcripts AND ERCC spike-ins, and
#empty sections at the beginning and the end of the dataset

ggplot(bc.stats.2.tomo13, aes(x = Section, y = Transcripts)) + 
  theme_bw() +
  geom_bar(stat="identity") +
  #geom_line(aes(x = Section, y = Transcripts, group=1)) +
  theme(axis.text.x = element_text(angle = 60, size=15),
        text=element_text(size=24),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Animal-to-vegetal position") + 
  ylab("Sequencing reads") +
  labs(title = "Total read distribution in on-cell stage embryo")
#ggsave("tomo13.read.profile.pdf", width = 16, height = 9, unit="in")

ggplot(bc.stats.2.tomo13, aes(x = Section, y = Genes, group=1)) + geom_line() +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45)) +
  labs(title = "Genes in sections with more than 8000 transcripts in tomo13")

gene.trans.tomo13 <- merge(add.names, transcripts.tomo13) #that converts GENEIDs into gene names
gene.trans.2.tomo13 <- gene.trans.tomo13[ ,-1]

length(unique(gene.trans.2.tomo13$Gene))

###Normalisation of data
# Filter for good sections, filter out lowly expressed genes and normalize expression to reads per sections for tomo13####
#filter out empty or low quality sections
transcripts.fs.tomo13 <- gene.trans.2.tomo13[, c(T, colSums(gene.trans.2.tomo13[, -1]) > min.transcripts)]
#filter out very lowly expressed genes - here, we demand 5 reads (min.e) in at least one section (min.s)
transcripts.ffs.tomo13 <- transcripts.fs.tomo13
transcripts.ffs.tomo13$o.max <- apply(transcripts.fs.tomo13[, -1], 1, function(x) sum(x > min.e))
transcripts.ffs.tomo13 <- transcripts.ffs.tomo13[transcripts.ffs.tomo13$o.max >= min.s, -ncol(transcripts.ffs.tomo13)]
# normalize to reads per section - divide every count by the reads in that section and multiply with a constant value (here the
# median of transcript recovery, Seurat uses 10^4)
transcripts.ffs.n.tomo13 <- transcripts.ffs.tomo13
median.tomo13 <- median(colSums(transcripts.ffs.n.tomo13[, -1]))
transcripts.ffs.n.tomo13[, -1] <- median.tomo13 * t(t(transcripts.ffs.n.tomo13[, -1])/colSums(transcripts.ffs.n.tomo13[, -1]))

colSums(transcripts.ffs.n.tomo13[,-1]) # -> constant

# Checking known localized mRNAs ####
#animally localized
plot.gene.name("exd2", transcripts.ffs.n.tomo13)
#ggsave("tomo13.exd2.pdf", dpi=300, width=7, height = 5, units = "in")

#vegetally located genes, this chunks produces plots as shown in figure 1c:
plot.gene("dazl", transcripts.ffs.n.tomo13)
plot.gene("dazl", gene.trans.2.tomo13)
plot.gene("celf1", transcripts.ffs.n.tomo13)
plot.gene("trim36", transcripts.ffs.n.tomo13)
plot.gene("wnt8a", transcripts.ffs.n.tomo13)
plot.gene("grip2a", transcripts.ffs.n.tomo13)

# Z-Score transition of normalized counts####
transcripts.z.tomo13 <- transcripts.ffs.n.tomo13
transcripts.z.tomo13[, -1] <- t(scale(t(transcripts.ffs.n.tomo13[, -1])))

plot.gene.name("vasa", transcripts.ffs.n.tomo13)
plot.gene.name("dazl", transcripts.ffs.n.tomo13)
plot.gene.name("grip2a", transcripts.ffs.n.tomo13)

plot.gene.name("vasa", transcripts.z.tomo13)
plot.gene.name("dazl", transcripts.z.tomo13)
plot.gene.name("grip2a", transcripts.z.tomo13)
#ggsave("tomo13.grip2a.z.pdf", dpi=300, width=7, height = 5, units = "in")

#Self-organizing maps for clustering of the genes, coded by Bastiaan Spanjaard####
# Calculate the cumulative expression for all genes, normalized to the maximum expression. Then cluster genes with SOM.
#colnames(transcripts.ffs.n.tomo13)=c("Gene", paste("X", 1:69, sep = "."))
transcripts.tomo13.cum <- plyr::adply(transcripts.ffs.n.tomo13, 1, 
                               function(x){
                                 y <- x
                                 y[-1] <- cumsum(t(x[-1]))/sum(x[-1])
                                 y <- data.frame(y)
                                 return(y)})

plot.gene.name("dazl", transcripts.tomo13.cum)
plot.gene.name("bmp1a", transcripts.tomo13.cum)
plot.gene.name("fanci", transcripts.tomo13.cum)
#ggsave("tomo13.bmp1a.cum.pdf", dpi=300, width=7, height = 5, units = "in")
#reproduces plots in supplemental Fig. 2b

# SOM of cumulative expression ####
# Calculate and plot the Self-Organizing Map in a 25x16-grid.
tomo13.som.grid <-
  supersom(data = as.matrix(transcripts.tomo13.cum[, -1]), grid = somgrid(25, 16, "rectangular"))

plot(tomo13.som.grid, "codes")

tomo13.som.grid <-
  supersom(data = as.matrix(transcripts.tomo13.cum[, -1]), grid = somgrid(50, 1, "rectangular"))

# Calculate and plot the Self-Organizing Map in a column
tomo13.som <- supersom(data = as.matrix(transcripts.tomo13.cum[, -1]), 
                  grid = somgrid(1, 50, "rectangular"))
tomo13.codes.melt <- melt(tomo13.som$codes, 
                          varnames = c("Profile", "Section"),
                          value.name = "Cum.expression")
tomo13.codes.melt$Section <- factor(tomo13.codes.melt$Section,
                                    levels = paste("X", 1:96, sep = ""))
tomo13.codes.melt$Section.num <- as.integer(sub("X", "", tomo13.codes.melt$Section))

#the following code chuck reproduces the plot in supplemental Figure 2c:
ggplot(tomo13.codes.melt) +
  geom_tile(aes(x = Section, y = Profile, fill = Cum.expression)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0.5) +
  scale_x_discrete(breaks = paste("X.", 25 + 10*(0:6), sep = ""),
                   labels = 25 + 10*(0:6)) +
  theme(panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Cumulative expression")+
  labs(title="Self-organizing map of cumulated gene expression")
#ggsave("tomo13_1_50_density_genes_cumulative_som.png", dpi=300, width=12, height = 7, units = "in")
#the SOM algorithm will randomly put vegetally localised genes in profiles 1-3 or 47-50, pay attention to that when using
#the next lines

# Determine and plot how many genes are in each profile
tomo13.codes.density <- data.frame(table(tomo13.som$unit.classif))

ggplot(tomo13.codes.density) +
  geom_tile(aes(x = 1, y = Var1, fill = Freq)) +
  scale_y_discrete(breaks = 20 * 1:5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 50) +
  labs(x = "", y = "", fill = "Frequency") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank())+
  labs(title="Gene density per profile")
#ggsave("tomo13_SOM50_density_per_profile.png", dpi=300, width=7, height = 5, units = "in")

# Assign which gene belongs to which profile
tomo13.som.clustering <- 
  data.frame(Gene = transcripts.tomo13.cum$Gene,
             Profile = tomo13.som$unit.classif,
             Transcripts = rowSums(transcripts.ffs.n.tomo13[, -1]))
tomo13.som.profile=tomo13.som.clustering[order(tomo13.som.clustering$Profile, decreasing = F),]

SOM.candidates=merge(tomo13.som.profile, add.names)
SOM.candidates=SOM.candidates[order(SOM.candidates$Profile, decreasing = F),]
#write.csv(SOM.candidates[SOM.candidates$Profile<48,], "tomo13.SOM50.neglist.csv", quote=F, row.names = F)
can=SOM.candidates[SOM.candidates$Profile>48,]
length(can$Gene)

#Plotting z-scores according to SOM profiles####
tomo13.z.SOM=merge(transcripts.z.tomo13, tomo13.som.clustering)
tomo13.z.SOM=tomo13.z.SOM[,-ncol(tomo13.z.SOM)]
tomo13.z.SOM=tomo13.z.SOM[,-(ncol(tomo13.z.SOM)-1)]
#sort genes according to SOM profile number
tomo13.z.SOM=tomo13.z.SOM[order(tomo13.z.SOM$Profile, decreasing = T),]

tomo13.z.SOM.melt= melt(tomo13.z.SOM, 
                          id.vars = c("Profile", "Gene"),
                          variable.name = c("Section.Nr"),
                          value.name = "Z.score")
tomo13.z.SOM.melt$Section = (sub("X", "", tomo13.z.SOM.melt$Section.Nr))

tomo13.z.SOM.melt.order=tomo13.z.SOM.melt[order(tomo13.z.SOM.melt$Profile, decreasing = T),]

tomo13.z.SOM.melt$Gene =
  factor(tomo13.z.SOM.melt$Gene, 
         levels = tomo13.z.SOM.melt$Gene[order(-tomo13.z.SOM.melt$Profile)])

require(scales)

#this code chuck reproduces the plot in Figure 2a:
ggplot(tomo13.z.SOM.melt)+
  geom_tile(aes(x = Section, y = Gene, fill = Z.score), show.legend=NA) +
    #scale_fill_gradient2(low = "steelblue4", high = "darkorange2", mid = "white") +
  scale_fill_gradientn(colours = c("steelblue4", "white", "darkorange2"),
                       limits = c(-2, 2)) +
  theme(text=element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Animal-to-vegetal position") + 
  labs(title="Transcript localization in a one cell stage embryo")
#ggsave("tomo13.SOM50.heatmap.orange.blue.png", width = 3.5, height = 2, unit="in")

#compares the gene recovery between different tomo-seq datasets:
detected.genes.tomo17 <- data.frame("Genes"=transcripts.tomo17[,1],
                                 "rep2"=rowSums(transcripts.tomo17[,-1]))

detected.genes.tomo13 <- data.frame("Genes"=transcripts.tomo13[,1],
                                 "rep1"=rowSums(transcripts.tomo13[,-1]))

#build a dataframe for all your replicates:
expr <- merge(detected.genes.tomo17, detected.genes.tomo13, by="Genes", all = T)

length(unique(expr$Genes))

sum(is.na(expr))
expr[is.na(expr)]=0.1 #set the NA values to a very small, arbitrary value -> find these values back at the axes, respectively
#normalize to sequencing depth
s <- as.numeric(colSums(expr[,c(2,3)]))
expr$rep2.norm <- expr$Tomo17/s[1]*10^6
expr$rep1.norm <- expr$Tomo13/s[2]*10^6
colSums(expr[,-1])

#calculate expression correlation
cor(expr$rep1.norm, expr$rep2.norm, method= "pearson")
#pearson= 0.993989

#this piece of code produces the figure 1e:
ggplot(data=expr,aes(x=rep2.norm, y=rep1.norm))+
  geom_point() +
  theme_bw() +
  scale_y_log10() + scale_x_log10() +
  annotate("text", label="pearsons R = 0.99", x=30, y=100000)+
  labs(title="Correlation of tomo-seq replicates")+
  #geom_text(aes(label=Genes), hjust="left", size= 3, check_overlap = T) +
  geom_smooth(method="lm", colour="darkorange")+
  xlab("replicate 2")+ ylab("replicate 1")

# Generate a list of vegetally localized genes  in three replicate samples:
high.conf=tomo.replicates.SOM[tomo.replicates.SOM$SOM.rep1>47&(tomo.replicates.SOM$SOM.rep2>45&tomo.replicates.SOM$SOM.rep2>45)|tomo.replicates.SOM$SOM.rep3>47&(tomo.replicates.SOM$SOM.rep1>45&tomo.replicates.SOM$SOM.rep2>45)|tomo.replicates.SOM$SOM.rep2>47&(tomo.replicates.SOM$SOM.rep1>45&tomo.replicates.SOM$SOM.rep3>45),]
high.conf=high.conf[order(high.conf$SOM.rep1, high.conf$SOM.rep3, high.conf$SOM.rep2,decreasing = T),]
#this results in a list of 97 confidentially vegetally localizing genes, find it as "vegetal.genes.dR.csv" in the tomo-seq repository on github

#Venn diagram for vegetally localized genes between replicates
#this pieces produces the plot in Figure 2b:
grid.newpage()
draw.triple.venn(area1 = length(unique(tomo.replicates.SOM$gene[tomo.replicates.SOM$SOM.rep3>47])),
                 area2 = length(unique(tomo.replicates.SOM$gene[tomo.replicates.SOM$SOM.rep1>47])),
                 area3 = length(unique(tomo.replicates.SOM$gene[tomo.replicates.SOM$SOM.rep2>47])),
                 n12 = 84, n23 = 71, n13 = 83, n123 = length(high.conf.consistent$Gene), category = c("rep3", "rep1", "rep2"), lty = "blank",
                 fill = c("turquoise", "yellow", "orange"))

#the following code generates the plot in 2c, correlating gene expression of genes in profiles 46 and higher (vegetally localized genes):
tomo.replicates.SOM[is.na(tomo.replicates.SOM)] = 0.1

cor.plot <-
  ggplot()+
  theme_bw() + 
  geom_point(data = tomo.replicates.SOM[tomo.replicates.SOM$SOM.rep2>45 | tomo.replicates.SOM$SOM.rep3 > 45,],
             aes(x= expr.rep2, y=expr.rep3),
             alpha = 0.5, size = 2) +
  scale_x_log10() + scale_y_log10() 

xplot <- ggplot()+
  theme_classic() +
  geom_density(data = tomo.replicates.SOM[tomo.replicates.SOM$SOM.rep2>45 | tomo.replicates.SOM$SOM.rep3 > 45,], 
               aes(expr.rep2), fill='lightgrey') + scale_x_log10() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
yplot <- ggplot() +
  theme_classic() +
  geom_density(data = tomo.replicates.SOM[tomo.replicates.SOM$SOM.rep2>45 | tomo.replicates.SOM$SOM.rep3 > 45,], 
               aes(expr.rep3), fill = "lightgrey")+ scale_x_log10() + 
  theme(axis.title.y =element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_flip()

#instead of ggarrange, use grid.arrange to align the marginal density plots. For that, create 
#a blank plot as a spaceholder

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

grid.arrange(xplot, blankPlot, cor.plot, yplot, 
             ncol = 2, nrow = 2, 
             widths = c(2, 0.5), heights = c(0.5, 2))
