library(tidyverse)
library(reshape2)
library(broom)
library(multcomp)
library(car)
library(ggrepel)
library(data.table)
library(phylosignal)
## Phylogenetic signal analysis

setwd("")
tree <- read.tree("branch_len.tree") # need branch lengths for phylogenetic signal detection
tree$tip.label
tree$tip.label <- c("chlmay", "ppatens", "smoel", "atricho", "osat", "zmays", "athal", 
                    "esals", "slyco", "stuber")
tree <- as(tree, "phylo4")


## import all lncrna info csvs!

setwd("~/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/script/lncrnainfo/")
temp = list.files(pattern="*lncinfo.csv")
lncinfo = lapply(temp, read.csv)
lncinfo <- lapply(lncinfo, function(x) x[,-1])
lncinfo <- lapply(lncinfo, "[", TRUE, -c(1))


speciesnames <- gsub("*_lncinfo.csv$", "", temp)
lncinfo <- setNames(lncinfo, speciesnames)



x <- melt(lncinfo, id=c("qry_id", "prediction", "class_code_x"))

x$value <- as.numeric(x$value)
x$prediction <- as.factor(x$prediction)
x$L1 <- as.factor(x$L1)
x$prediction <- as.factor(x$prediction)

all_features <- c("num_exons_y", "ORF", "len_y", "GC")

x_sub <- x %>% filter(., variable %in% all_features)



# summarizing data to add to phylogeny
phylo_data <- x %>%
  group_by(L1, variable, prediction) %>%
  summarize(average = mean(value, na.rm = TRUE))
filt_data <- phylo_data %>% filter(variable %in% all_features)
lnc_data <- filt_data %>% filter(prediction == 1)
nonlnc_data <- filt_data %>% filter(prediction ==0)



lnc_tree <- dcast(lnc_data, L1 ~ variable, value.var = 'average')

nonlnc_tree <- dcast(nonlnc_data, L1 ~ variable, value.var = 'average')

species_order <- c("chlamy", "ppatens", "smoel", "atricho", "osat", "zmays", "athal", "esals", "slyco", "stuber")


lnc_tree <- lnc_tree[match(species_order, lnc_tree$L1),]
nonlnc_tree <- nonlnc_tree[match(species_order, nonlnc_tree$L1),]


rownames(lnc_tree) <- NULL
rownames(nonlnc_tree) <- NULL

diff <- nonlnc_tree[,-1]-lnc_tree[,-1]


## Getting data into phylo4d format
tree4d_lnc <- phylo4d(tree, tip.data=lnc_tree[,-1])
tree4d_nonlnc <- phylo4d(tree, tip.data=nonlnc_tree[,-1])
tree4d_diff <- phylo4d(tree, tip.data=diff)
barplot(tree4d_diff)
barplot(tree4d_nonlnc)


## Getting signal estimates
set.seed(2248)
diff_signal <- phyloSignal(tree4d_diff, methods=c("Lambda", "K", "I"), reps=1000)
round(diff_signal$stat,digits=3 )## REQUIRES BRANCH LENGTHS!!
round(diff_signal$pvalue, digits = 4)
set.seed(2248)
lnc_signal <- phyloSignal(tree4d_lnc,  methods=c("Lambda", "K", "I"), reps=1000)
round(lnc_signal$stat, digits=3)
round(lnc_signal$pvalue, digits = 4)

## locating where in tree there is signal
set.seed(2248)
orf.corr <- phyloCorrelogram(tree4d_lnc, trait = "ORF")
plot.correlogram(orf.corr)
gc.corr <- phyloCorrelogram(tree4d_lnc, trait = "GC")
length.corr <- phyloCorrelogram(tree4d_lnc, trait = "len_y")
exon.corr <- phyloCorrelogram(tree4d_lnc, trait="num_exons_y")
lnc.corr.all <- phyloCorrelogram(tree4d_lnc, trait=c("ORF", "GC", "len_y", "num_exons_y"))
lnc.corr.three <- phyloCorrelogram(tree4d_lnc, trait=c("ORF", "GC","num_exons_y"))
plot(lnc.corr)


## extracting data from phyloCorrelogram for custom plotting
correlogram_cust <- function(gc.corr, orf.corr, length.corr, exon.corr){
  h0 <- -1/(10-1) #h0 is the value we expect without phylogenetic signal
  gc<- data.frame(gc.corr$res, trait=rep("GC", 100)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Not Significant"))
  orf <- data.frame(orf.corr$res, trait=rep("ORF", 100))  %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Not Significant"))
  length <- data.frame(length.corr$res, trait=rep("length", 100)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait"))  %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Not Significant"))
  exon <- data.frame(exon.corr$res, trait=rep("exon", 100)) %>% 
    `colnames<-`(c("dist", "min", "max", "cor", "trait")) %>% 
    mutate(sig = case_when(min > h0 & max > h0 ~ "Significant",
                           min < h0 & max < h0 ~ "Significant", 
                           TRUE ~ "Not Significant"))
  plotdata <- rbind(gc, orf, length, exon)
  return(plotdata)
  }

lnccor.plotdata <- correlogram_cust(gc.corr=gc.corr, orf.corr = orf.corr, length.corr= length.corr, 
                         exon.corr=exon.corr)

ggplot(lnccor.plotdata, aes(x = dist, y = cor, color = sig, group = trait)) + 
  geom_point(aes(shape=sig)) + geom_line() + 
  #geom_line(aes(y=min), color = "black", alpha=0.3) +
  #geom_line(aes(y=max), color="black", alpha=0.3) +
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3, color = NA) +
  facet_grid(~trait) +
  geom_hline(yintercept=(-1/9), color="black") 
 
lnccor.plotdata <- correlogram_cust(gc.corr=gc.corr, orf.corr = orf.corr, length.corr= length.corr, 
                                    exon.corr=exon.corr)



lnc.lipa <- lipaMoran(tree4d_lnc)
lnc.lipa.p4d <- lipaMoran(tree4d_lnc, as.p4d = TRUE)
barplot.phylo4d(tree4d_lnc, bar.col=(lnc.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE)


set.seed(2248)
nonlnc_signal <- phyloSignal(tree4d_nonlnc, methods=c("Lambda", "K", "I"), reps=1000)
round(nonlnc_signal$stat, digits=3)
round(nonlnc_signal$pvalue, digits = 4)

orf.corr.non <- phyloCorrelogram(tree4d_nonlnc, trait = "ORF")
gc.corr.non <- phyloCorrelogram(tree4d_nonlnc, trait = "GC")
length.corr.non <- phyloCorrelogram(tree4d_nonlnc, trait = "len_y")
exon.corr.non <- phyloCorrelogram(tree4d_nonlnc, trait="num_exons_y")
nonlnc.corr.all <- phyloCorrelogram(tree4d_nonlnc, trait=c("ORF", "GC", "len_y", "num_exons_y"))
nonlnc.corr.three <- phyloCorrelogram(tree4d_nonlnc, trait=c("ORF", "GC","num_exons_y"))

noncor.plotdata <- correlogram_cust(gc.corr=gc.corr.non, orf.corr = orf.corr.non, 
                                    length.corr= length.corr.non,exon.corr=exon.corr.non)


ggplot(noncor.plotdata, aes(x = dist, y = cor, color = sig, group = trait)) + 
  geom_point(aes(shape=sig)) + #geom_line() + 
  #geom_line(aes(y=min), color = "black", alpha=0.3) +
  #geom_line(aes(y=max), color="black", alpha=0.3) +
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.3, color = NA) +
  facet_grid(~trait) +
  geom_hline(yintercept=(-1/9), color="black") 

all.plotdata <- rbind(data.frame(noncor.plotdata, type=rep("nonlnc", 400)), 
                           data.frame(lnccor.plotdata, type=rep("lnc", 400)))

all.plotdata$trait <- factor(all.plotdata$trait, levels = c("ORF", "GC", "exon", "length"))
all.plotdata$sig <- factor(all.plotdata$sig, levels = c("Significant", "Not Significant"))
legend_ord <- levels(with(all.plotdata, reorder(c("Significant", "Not Significant"))))
                           
feature_names <- as_labeller(c(
  'GC'="GC%",
  'length'="Transcript length (nt)",
  'exon'="Number of exons",
  'ORF'="ORF length (nt)",
  'lnc'="lncRNA",
  'nonlnc' = "All other transcripts"
))


correlogram <- ggplot(all.plotdata, aes(x = dist, y = cor, color = sig, group = trait)) + 
  geom_point() + #geom_line() + 
  #geom_line(aes(y=min), color = "black", alpha=0.3) +
  #geom_line(aes(y=max), color="black", alpha=0.3) +
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.15, color = NA) +
  facet_grid(type ~ trait, labeller=feature_names) +
  geom_hline(yintercept=(-1/9), color="black") + theme_bw() + ylab("Correlation") +
  xlab("Phylogenetic distance") + scale_color_manual(values=c("#ed7a53","#A7c6da")) +#,
                                                    # guide = guide_legend(reverse = TRUE)) + 
  theme(legend.title=element_blank(), legend.position = "bottom")

ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/figs/correlogram.pdf",
       correlogram)  
