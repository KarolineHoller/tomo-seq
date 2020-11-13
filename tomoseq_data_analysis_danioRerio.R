# Description ####
# How we analyzed tomoSeq data for zebrafish 1-cell stage embryos (here:tomo13)

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

setwd("~/data/junker/users/Karo/collab_Meyer/count_tables/danio_tomoseq_raw_counts")
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


# Functions ####
find.peak <- function(gene.expression, peak.length, cutoff){
  s = 0
  peak.found = F
  for(i in 1:length(gene.expression)){
    if(gene.expression[i] > cutoff){
      s = s + 1
      if(s >= peak.length){
        peak.found <- T
      }
    }else{
      s <- 0
    }
  }
  return(peak.found)
}

overseq <- function(x, y, name){
  os <- as.matrix((x+0.01)/(y+0.01))
  #os[is.infinite(os)] <- 0
  #os[is.na(os)] <- 0
  os.melt <- melt(os)
  print(
    ggplot(os.melt, aes(x = value)) + theme_light() +
      geom_histogram(binwidth = 0.2) +
      scale_y_log10() +
      labs(title = paste("sequencing saturation", name),
           x = "read counts / UMI counts") #+
      #scale_x_continuous(breaks = 1:35)
  )  
}

plot.gene.line <- function(gene.expression){
  plot.title <- paste("Spatial expression of", gene.expression[1])
  plot.data <- 
    data.frame(Section = as.integer(sapply(colnames(gene.expression)[-1], function(x) unlist(strsplit(x, "X"))[2])),
               Expression = t(gene.expression[-1]))
  colnames(plot.data)[2] <- "Expression"
  print(
    ggplot(plot.data, aes(x = Section, y = Expression)) + geom_line() +
      labs(title = plot.title)
  )
}

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

plot.short <- function(name, expression.data){
  gene.expression <- expression.data[expression.data$Gene == name, ]
  plot.data <- 
    data.frame(Section = as.integer(sapply(colnames(gene.expression)[-1], function(x) unlist(strsplit(x, "X"))[2])),
               Expression = t(gene.expression[-1]))
  colnames(plot.data)[2] <- "Expression"
  p=ggplot(plot.data, aes(x = Section, y = Expression)) +geom_bar(stat = "identity") +
    theme_bw() +
    theme(text = element_text(size=6),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank())+
    xlab("")+
    ylab("")
  return(p)
}

# Load data ####
#the here mentioned files are concatenated count tables from 2 consequtive sequencing runs. The
# librariea are from a tomo-seq experiment from 1-cell stage zebrafish embyros

transcripts.in.tomo13=read.table("~/data/junker/users/Karo/collab_Meyer/count_tables/danio_tomoseq_raw_counts/tomo13_deep.coutt.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(transcripts.in.tomo13)[-1]=paste("X", 1:96, sep = "")
transcripts.tomo13 = transcripts.in.tomo13[!grepl("ERCC", transcripts.in.tomo13$GENEID), ]
counts.in.tomo13 <- read.table("~/data/junker/users/Karo/collab_Meyer/count_tables/danio_tomoseq_raw_counts/tomo13_deep.coutc.csv", sep = "\t", header = T, stringsAsFactors = F)
counts.tomo13 <- counts.in.tomo13[!grepl("ERCC", counts.in.tomo13$GENEID), ]
UMIs.in.tomo13 <- read.table("~/data/junker/users/Karo/collab_Meyer/count_tables/danio_tomoseq_raw_counts/tomo13_deep.coutb.csv", sep = "\t", header = T, stringsAsFactors = F)
UMIs.tomo13 <- UMIs.in.tomo13[!grepl("ERCC", UMIs.in.tomo13$GENEID), ]
add.names= read.table("~/data/junker/users/Karo/collab_Meyer/count_tables/danio_tomoseq_raw_counts/gene_names_translator.csv", sep = ",", header=T, stringsAsFactors = F)
#downloaded from ensemble, matching the transcriptome that we mapped to

#ERCC content per section in tomo13 - checks if the internal controls are distributed evenly (few dropouts are marks of high quality)
ERCC.tomo13=transcripts.in.tomo13[grep("ERCC", transcripts.in.tomo13$GENEID),]
#colnames(ERCC.tomo13)
ERCC.tomo13.s=colSums(ERCC.tomo13[,-1])
plot(ERCC.tomo13.s, type="h",
     xlab="Section number", ylab="ERCC reads",
     main="ERCC reads in tomo13")

# Sequencing statistics tomo13####

#the number of reads per UMI are a measure of how deeply a library is sequenced - desirable would be that every UMI has been seen multiple times
overseq(counts.tomo13[, -1], UMIs.tomo13[, -1], "tomo-seq")

df_overseq <- data.frame("gene" = counts.tomo13$GENEID,
                         "reads" = rowSums(counts.tomo13[,-1]),
                         "UMI" = rowSums(UMIs.tomo13[,-1]))

#function for UMI saturation (Grün et al 2014):
#while 4^6*96 is the maximum amount of UMIS that we could potentially recover, with 6 bp random
#barcode, multiplied with 96 sections
UMI.saturation <- function(x) {
  (1- exp(-x/(4^6*96)))*4^6*96
}

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

#visualizes the distribution of transcript- and read recovery per section
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
  labs(title = "Genes in sections with more than 1000 transcripts in tomo13")
#ggsave("genes_per_f_section_tomo13.cut1000.pdf", width = 8, height = 5, unit="in")

gene.trans.tomo13 = merge(add.names, transcripts.tomo13) #that converts GENEIDs into gene names
gene.trans.2.tomo13 = gene.trans.tomo13[ ,-1]

length(unique(gene.trans.2.tomo13$Gene))
#16879
#write.csv(gene.trans.3.tomo13, "tomo13_wo_filtering.csv",
 #         quote = F, row.names = F)

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
#write.csv(transcripts.ffs.n.tomo13, "tomo13_norm2section_transcripts.csv",
#          quote = F, row.names = F)

# Checking known localized mRNAs ####
#animally distributed
plot.gene.name("pou2f1b", transcripts.ffs.n.tomo13)
plot.gene.name("bmp1a", transcripts.ffs.n.tomo13)
plot.gene.name("pabpn1l", transcripts.ffs.n.tomo13)
plot.gene.name("sox19", transcripts.ffs.n.tomo13)
plot.gene.name("CKAP2", transcripts.ffs.n.tomo13)
plot.gene.name("cth1", transcripts.ffs.n.tomo13)
plot.gene.name("exd2", transcripts.ffs.n.tomo13)
plot.gene.name("nkap", transcripts.ffs.n.tomo13)
#ggsave("tomo13.nkap.pdf", dpi=300, width=7, height = 5, units = "in")

#vegetally located genes
plot.gene.name("dazl", transcripts.ffs.n.tomo13)
plot.gene.name("dazl", gene.trans.2.tomo13)
plot.gene.name("celf1", transcripts.ffs.n.tomo13)
plot.gene.name("wnt8a", transcripts.ffs.n.tomo13)
plot.gene.name("grip2a", transcripts.ffs.n.tomo13)

# Z-Score transition ####
transcripts.z.tomo13 = transcripts.ffs.n.tomo13
transcripts.z.tomo13[, -1] <- t(scale(t(transcripts.ffs.n.tomo13[, -1])))

plot.gene.name("vasa", transcripts.ffs.n.tomo13)
plot.gene.name("dazl", transcripts.ffs.n.tomo13)
plot.gene.name("grip2a", transcripts.ffs.n.tomo13)

plot.gene.name("vasa", transcripts.z.tomo13)
plot.gene.name("dazl", transcripts.z.tomo13)
plot.gene.name("grip2a", transcripts.z.tomo13)
#ggsave("tomo13.grip2a.z.pdf", dpi=300, width=7, height = 5, units = "in")

#Self-organizing maps for clustering of the genes, Bastiaan Spanjaard####
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


# SOM of cumulative expression ####
# Calculate and plot the Self-Organizing Map in a 25x16-grid.
tomo13.som.grid <-
  supersom(data = as.matrix(transcripts.tomo13.cum[, -1]), grid = somgrid(25, 16, "rectangular"))

plot(tomo13.som.grid, "codes")
#ggsave("tomo13_gene_profiles_som.png", dpi=300, width=16, height = 9, units = "in")

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
tomo13.codes.melt$Section.num = as.integer(sub("X", "", tomo13.codes.melt$Section))

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

#For visualization, plot a high number of genes from the SOM list ####
#Plot the localization of the highest ranked vegetally localised 50 genes from the SOM list
p1=plot.short("celf1", transcripts.ffs.n.tomo13)
p2=plot.short("slc38a7", transcripts.ffs.n.tomo13)
p3=plot.short("dazl", transcripts.ffs.n.tomo13)
p4=plot.short("fgd6", transcripts.ffs.n.tomo13)
p5=plot.short("shcbp1", transcripts.ffs.n.tomo13)
p6=plot.short("sybu", transcripts.ffs.n.tomo13)
p7=plot.short("trim36", transcripts.ffs.n.tomo13)
p8=plot.short("cdk6", transcripts.ffs.n.tomo13)
p9=plot.short("camk2g1", transcripts.ffs.n.tomo13)
p10=plot.short("grip2a", transcripts.ffs.n.tomo13)
p12=plot.short("tpd52l2b", transcripts.ffs.n.tomo13)
p13=plot.short("sh2d5", transcripts.ffs.n.tomo13)
p14=plot.short("myl12.1", transcripts.ffs.n.tomo13)
p15=plot.short("hmmr", transcripts.ffs.n.tomo13)
p16=plot.short("clic4", transcripts.ffs.n.tomo13)
p17=plot.short("lmbr1", transcripts.ffs.n.tomo13)
p18=plot.short("UBL7", transcripts.ffs.n.tomo13)
p19=plot.short("fbxw11a", transcripts.ffs.n.tomo13)
p20=plot.short("lef1", transcripts.ffs.n.tomo13)
p21=plot.short("ndel1b", transcripts.ffs.n.tomo13)
p22=plot.short("sulf1", transcripts.ffs.n.tomo13)
p23=plot.short("GIT1", transcripts.ffs.n.tomo13)
p24=plot.short("gnav1", transcripts.ffs.n.tomo13)
p25=plot.short("KIF2A_$9297", transcripts.ffs.n.tomo13)
p26=plot.short("ptprga", transcripts.ffs.n.tomo13)
p27=plot.short("ca7", transcripts.ffs.n.tomo13)
p28=plot.short("pawr", transcripts.ffs.n.tomo13)
p29=plot.short("wnt8a_$10565", transcripts.ffs.n.tomo13)
p30=plot.short("CCDC88C", transcripts.ffs.n.tomo13)
p31=plot.short("slc7a6", transcripts.ffs.n.tomo13)
p32=plot.short("rftn2", transcripts.ffs.n.tomo13)
p33=plot.short("rab3aa", transcripts.ffs.n.tomo13)
p34=plot.short("mpzl1l", transcripts.ffs.n.tomo13)
p35=plot.short("vrtn", transcripts.ffs.n.tomo13)
p36=plot.short("ctdsplb", transcripts.fs.n.tomo13)
p37=plot.short("ino80db", transcripts.ffs.n.tomo13)
p38=plot.short("golga3", transcripts.ffs.n.tomo13)
p39=plot.short("zgc:158856", transcripts.ffs.n.tomo13)
p40=plot.short("fndc3ba", transcripts.ffs.n.tomo13)
p41=plot.short("ralbp1", transcripts.ffs.n.tomo13)
p42=plot.short("CU467832.1", transcripts.ffs.n.tomo13)
p43=plot.short("HOOK2", transcripts.ffs.n.tomo13)
p44=plot.short("si:dkey-121h17.7", transcripts.ffs.n.tomo13)

grid.arrange(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,
             p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36, ncol=7, nrow =5)
#ggsave("tomo13_SOM35_distr.bw.2.pdf", width=9, height=6, unit="in")
