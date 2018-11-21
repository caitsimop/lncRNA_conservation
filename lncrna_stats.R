## lncrna stats
library(tidyverse)
library(reshape2)
library(broom)
library(multcomp)
library(car)
library(ggrepel)
library(reshape2)
library(patchwork)
library(ggtree)
library(ape)

## import all lncrna info csvs!

setwd("~/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/script/lncrnainfo/")  ## set directory containing CSVs found in  https://github.com/caitsimop/lncRNA-compGenomics
temp = list.files(pattern="*lncinfo.csv")
lncinfo = lapply(temp, read.csv)
lncinfo <- lapply(lncinfo, function(x) x[,-1])
lncinfo <- lapply(lncinfo, "[", TRUE, -c(1))

speciesnames <- gsub("*_lncinfo.csv$", "", temp)
lncinfo <- setNames(lncinfo, speciesnames)


# melting data into tidyverse format
x <- melt(lncinfo, id=c("qry_id", "prediction", "class_code_x"))
x$value <- as.numeric(x$value)
x$prediction <- as.factor(x$prediction)



x$L1 <- as.factor(x$L1)
x$prediction <- as.factor(x$prediction)
### t-test to compare means within each species
 variable_ttest<-function(species_name,character, df){
 set.seed(123)
    df %>% 
    filter(variable == character & L1 == species_name) %>%
    t.test(value ~ prediction, data = .) %>% 
    tidy(.)%>% 
    add_column(., species = species_name) %>%
    add_column(., feature = character) %>%
    rbind(all_tests, .)
}
 
all_features <- c("num_exons_y", "ORF", "len_y", "GC")
all_tests <- data.frame(estimate = as.numeric(), estimate1= as.numeric(), estimate2 = as.numeric(),
                        statistic=as.numeric(),p.value=as.numeric(), parameter = as.numeric(), 
                        conf.low=as.numeric(), conf.high=as.numeric(),method=as.character(), 
                        alternative=as.character(), species=as.character(), feature=as.character())
for (n in speciesnames){
  for (i in all_features){
    all_tests <- variable_ttest(species_name=n, character=i, df=x)
  }}
i
## USING A NON-PARAMETRIC TEST!!
all_tests$adj.pval <- p.adjust(all_tests$p.value, method="fdr", n=nrow(all_tests))
all_tests[which(all_tests$adj.pval >= 0.05),]


## paired T tests of all species... 
filt_tests <- all_tests %>% filter(feature %in% all_features)

melted_x <- x

## let's summarize our actual data for nice boxplots!
summarized <- melted_x %>% filter(variable %in% all_features) %>%
  group_by(L1, prediction, variable) %>%
  summarize(
    count = n(),
    mean = mean(value),
    sd = sd(value),
    sem = sd(value)/sqrt(length(value)),
    ymin= min(value),
    ymax=max(value),
    upper=mean(value)+(sd(value)/sqrt(length(value))),
    lower=mean(value)-(sd(value)/sqrt(length(value)))) 

MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

feature_names <- list(
  'GC'="GC %",
  'len_y'="Transcript length",
  'num_exons_y'="Number of exons",
  'ORF'="ORF length"
)

## 2 way ANOVA!!
## lets do phylogenetic ANOVA
library(geiger)
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
rownames(diff) <- species_order

## ANOVA without accounting for 
ORF_anova <- x %>% filter(variable=="ORF") %>%
  aov(value ~ L1*prediction, data=.) 
drop1(ORF_anova, test="F")) ## smaller AIC without removing the interaction...thus the interaction should st
tidy(ORF_anova)

GC_anova <- x %>% filter(variable=="GC") %>%
  aov(value ~ L1*prediction, data=.) 
drop1(GC_anova, test="F") # keep interaction

exon_anova <- x %>% filter(variable=="num_exons_y") %>%
  aov(value ~ L1*prediction, data=.)
drop1(exon_anova, test="F") # keep interaction

len_anova <- x %>% filter(variable=="len_y") %>%
  aov(value ~ L1*prediction, data=.)
drop1(len_anova, test="F") # keep interaction

library(emmeans)


## ppost hoc tests to determine which species varied (ANOVA)
comp_ORF <- tidy(TukeyHSD(ORF_anova, "L1"))
comp_ORF[which(comp_ORF$adj.p.value > 0.05),]

#TukeyHSD(ORF_anova, "prediction")
comp_GC <- tidy(TukeyHSD(GC_anova, "L1"))
comp_GC[which(comp_GC$adj.p.value > 0.05),] #all different!!

comp_exon <- tidy(TukeyHSD(exon_anova, "L1"))
comp_exon[which(comp_exon$adj.p.value > 0.05),]

comp_len<- tidy(TukeyHSD(len_anova, "L1"))
comp_len[which(comp_len$adj.p.value > 0.05),]


## Figures for paper
all_tests$label <- c("", "a", "", "", 
                     "", "b", "a", "",
                    "", "", "", "",
                    "", "a", "b", "", 
                    "a", "c", "b", "", 
                    "", "", "", "", 
                    "b,c", "c", "a", "", 
                    "a,b", "", "", "",
                    "", "b", "", "", 
                    "c", "", "", "")

species_order <- c("chlamy", "ppatens", "smoel", "atricho", "osat", "zmays", "athal",
                   "esals", "slyco", "stuber")

## Change species nicknames to actual names
all_tests2 <- all_tests %>% mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                               osat = "O. sativa",
                               chlamy = "C. reinhardtii",
                               ppatens = "P. patens",
                               esals = "E. salsugineum",
                               slyco = "S. lycopersicum", 
                               athal = "A. thaliana", 
                               smoel = "S. moellendorffi",
                               zmays = "Z. mays",
                               stuber = "S. tuberosum"))

filt_data$L1 <- factor(filt_data$L1, levels = species_order) 
x2 <- filt_data %>% mutate(species = dplyr::recode(L1, atricho = "A. trichopoda",
                                                   osat = "O. sativa",
                                                   chlamy = "C. reinhardtii",
                                                   ppatens = "P. patens",
                                                   esals = "E. salsugineum",
                                                   slyco = "S. lycopersicum", 
                                                   athal = "A. thaliana", 
                                                   smoel = "S. moellendorffi",
                                                   zmays = "Z. mays",
                                                   stuber = "S. tuberosum"))  %>% 
                mutate(prediction = dplyr::recode(prediction, '1' = "lncRNA",
                                                 '0' = "All other transcripts"))


tree <- read.tree("/Users/caitlinsimopoulos/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/script/branch_len.tree")
tree$tip.label
tree$tip.label <- c("chlmay", "ppatens", "smoel", "atricho", "osat", "zmays", "athal", 
                    "esals", "slyco", "stuber")

feature_names <- as_labeller(c(
  'GC'="GC%",
  'len_y'="Transcript length (nt)",
  'num_exons_y'="Number of exons",
  'ORF'="ORF length (nt)"
))

independence <- c("#574D68")
sandstorm <- c("#F3D34A")
carmine <- c("#B0413E")

ggtreeplot <- ggtree(tree, branch.length="none")# + geom_tiplab()
lb <- get.tree(tree)$tip.label
d <- data.frame(label=lb, genus = c("Chlamydomomas", "Physcomitrella", 
                                     "Selaginella", "Amborella", "Oryza",
                                     "Zea", "Arabidopsis", "Eutrema", "Solanum",
                                     "Solanum"),
                species= c("reinhardtii", "patens", "moellendorffi","trichopoda", 
                           "sativa", "mays", "thaliana", "salsugineum", "lycopersicum",
                           "tuberosum"))

y <- ggtree(tree, branch.length="none") %<+% 
  d + geom_tiplab(aes(label=paste0('italic(', genus, ')~italic(', species, ')')), parse=TRUE) +
  xlim(0,15)
branch_len <- ggtree(tree) %<+% 
  d + geom_tiplab(aes(label=paste0('italic(', genus, ')~italic(', species, ')')), parse=TRUE) +
  xlim(0,15)

## The previous code plots mean trait values in the order of the phylogenetic tree, then combines the plots into one

ggsave("yourfilepath.pdf",
       plot=y, width = 9 , height = 4 , units = "in")

p <- ggplot(x2, aes(x=species, y=average, fill=prediction)) + 
  geom_point(stat = "identity",pch=21, size=3, colour="black") +
  coord_flip() +
  #geom_hline(yintercept=0, linetype="dashed", size=0.5) + guides(colour=FALSE) +
  facet_wrap(~variable, nrow=1, scales="free_x", labeller=feature_names) +
  labs(x = "", y = "Mean trait values") + scale_fill_manual(values = c("#D1CED5", sandstorm)) +
  theme(axis.text.y = element_text(face = "italic"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title=element_blank(),
        legend.position = "bottom")

p <- ggplot(subset(x2,prediction=="lncRNA"), aes(x=species, y=average)) + 
  geom_point(stat = "identity", pch = 21,size=4, colour="white", fill="#4E598C") +
  geom_point(data = subset(x2, prediction!="lncRNA"), aes(x=species, y=average, group = prediction), 
             fill=independence, pch=21, size=4, alpha=0.4, colour="white") +
  coord_flip() +
  #geom_hline(yintercept=0, linetype="dashed", size=0.5) + guides(colour=FALSE) +
  facet_wrap(~variable, nrow=1, scales="free_x", labeller=feature_names) + 
  #scale_fill_manual("predictions", values = c("lncRNA" = sandstorm, "All other transcripts" = independence)) +
  labs(x = "", y = "Mean trait values") + 
  #scale_fill_manual(name = "", values = c("point.label" = "#4E598C")) +
  #scale_color_manual(name = "", values = c("point.label" = "independence")) +
  scale_fill_identity(name = '', guide = 'legend',labels = c('lncRNA', 'All other transcripts')) +
  theme(axis.text.y = element_text(face = "italic"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title=element_blank(),
        legend.position = "bottom") #+ theme_light()


tree_plot <- ggtreeplot+ p + plot_layout(nrow = 1, widths = c(1, 4))
  
ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/figs/lnc_feat4.pdf",
       plot=tree_plot, width = 9 , height = 4 , units = "in")
## Different itertions of the figure.

### Barchart of novel # of lncrnas
lncperc <- c(5.14, 6.19, 2.98, 4.75, 7.55, 5.64, 16.55, 6.82, 5.24, 6.44, 3.07)
mean(lncperc)
std <- function(x) sd(x)/sqrt(length(x))
std(lncperc)

novellnc <- data.frame(read.csv("../novellncnum.csv", row.names = 1)) # a csv that contains numebr values for lncRNA counts
novellnc <- as.data.frame(t(novellnc))
novellnc$perc <- novellnc$novel/novellnc$Total
novellnc$species <- rownames(novellnc)
novellnc$oldperc <- novellnc$old/novellnc$Total
perc_novel <- novellnc %>% melt(., id=c("species", "Total")) %>% 
  filter(variable == "perc")
perc_novel$perc <- paste0(round(perc_novel$value*100, 1), "%") 

perc_novel <- perc_novel %>%
  mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                 osat = "O. sativa",
                                 chlamy = "C. reinhardtii",
                                 pppatens = "P. patens",
                                 esals = "E. salsugineum",
                                 slyco = "S. lycopersicum", 
                                 athal = "A. thaliana", 
                                 smoel = "S. moellendorffi",
                                 zmays = "Z. mays",
                                 stuber = "S. tuberosum"))


species_order <- c("stuber", "slyco", "esals", "athal", "zmays", "osat", "atricho", 
                   "smoel", "pppatens", "chlamy")
#species_order <- c("chlamy", "ppatens", "smoel", "atricho", "osat", "zmays", "athal",
#                   "esals", "slyco", "stuber")

all_tests$species <- factor(all_tests$species, levels = species_order)
levels(all_tests$species) <- c("S. tuberosum", "S. lycopersicum","E. salsugineum", "A. thaliana",
                               "Z. mays", "O. sativa", "A trichopoda", "S. moellendorffi", "P. patens",
                               "C. reinhardtii")

novel_bar_data <- novellnc %>% melt(., id=c("species")) %>% 
  filter(variable %in% c("old","novel")) %>%
  mutate(variable = dplyr::recode(variable, old = "In ref annotation",
                            novel = "Novel lncRNA"))
novel_bar_data$species <- factor(novel_bar_data$species, levels = species_order)


novel_bar <- novel_bar_data %>% mutate(variable = dplyr::recode(variable, old = "In ref annotation",
                                  novel = "Novel lncRNA")) %>%
  mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                 osat = "O. sativa",
                                 chlamy = "C. reinhardtii",
                                 pppatens = "P. patens",
                                 esals = "E. salsugineum",
                                 slyco = "S. lycopersicum", 
                                 athal = "A. thaliana", 
                                 smoel = "S. moellendorffi",
                                 zmays = "Z. mays",
                                 stuber = "S. tuberosum")) %>%
  ggplot(aes(y=value, x=species)) + 
  geom_bar(aes(fill=variable), color = "black", stat="identity") +
  theme_classic() + scale_fill_manual(values=c("white", "#ed7a53")) +
  geom_text(data = perc_novel, aes(label = perc,y = Total , x = species, vjust=-0.5)) +
  labs(x="Species", y="# of predicted lncRNAs") + 
  guides(fill=guide_legend(title="")) + theme(legend.position = c(0.2, 0.9), 
                                              axis.text.x = element_text(face = "italic", angle = 45, hjust = 1 ))
  
ggsave("~/Dropbox/McMaster/PhD/Written/Manuscripts/lncRNAapp_manuscript/figs/novel_bar.pdf",
       plot=novel_bar)

perc_old <- novellnc %>% melt(., id=c("species", "Total")) %>% 
  filter(variable == "oldperc")
perc_old$perc <- paste0(round(perc_old$value*100, 1), "%") 

old_bar_data <- novellnc %>% melt(., id=c("species")) %>% 
  filter(variable == "total")
old_bar_data$species <- factor(old_bar_data$species, levels = species_order)
perc_old <- perc_old %>%
  mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                 osat = "O. sativa",
                                 chlamy = "C. reinhardtii",
                                 pppatens = "P. patens",
                                 esals = "E. salsugineum",
                                 slyco = "S. lycopersicum", 
                                 athal = "A. thaliana", 
                                 smoel = "S. moellendorffi",
                                 zmays = "Z. mays",
                                 stuber = "S. tuberosum"))
perc_old$correct <- c("16.6%", "5.6%", "6.4%", "5.2%", "3.0%", "6.2%", "4.8%", "6.8%", "7.6%", "5.1%")

ref_info <- novellnc %>% melt(., id=c("species")) %>% 
  filter(variable %in% c("old","novel")) %>%
  mutate(variable = dplyr::recode(variable, old = "Predicted from ref. annotation",
                                  novel = "Predicted from RNASeq"))
ref_info$species <- factor(ref_info$species, levels = species_order)

ref_bar <- ref_info %>% mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                 osat = "O. sativa",
                                 chlamy = "C. reinhardtii",
                                 pppatens = "P. patens",
                                 esals = "E. salsugineum",
                                 slyco = "S. lycopersicum", 
                                 athal = "A. thaliana", 
                                 smoel = "S. moellendorffi",
                                 zmays = "Z. mays",
                                 stuber = "S. tuberosum")) %>%
  filter(variable == "Predicted from ref. annotation") %>%
  ggplot(aes(y=value, x=species)) + 
  geom_bar(aes(fill=variable), color = "black", stat="identity") +
  theme_classic() + #scale_fill_manual(values=c("#C9E4CA")) + #lavenderblush is ok too
  scale_fill_manual(values=c("#FDEF8D")) +
  #geom_text(data = perc_old, aes(label = correct,y = Total , x = species, vjust=-0.5), size=7) +
  labs(x="", y="# of predicted lncRNAs") + 
  guides(fill=guide_legend(title="")) + theme(legend.position = c(0.2, 0.9), 
                                              axis.text.x = element_text(face = "italic", 
                                              angle = 45, hjust = 1 ),text = element_text(size=30)) +
  scale_y_continuous(limits = c(0, 8000))
ggsave("~/Documents/jobapps/figeys/talk/figs/ref_bar.pdf",
       plot=ref_bar, width = 15, heigh = 8.6, units = "in")

perc_bar <- ref_info %>% mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                                        osat = "O. sativa",
                                                        chlamy = "C. reinhardtii",
                                                        pppatens = "P. patens",
                                                        esals = "E. salsugineum",
                                                        slyco = "S. lycopersicum", 
                                                        athal = "A. thaliana", 
                                                        smoel = "S. moellendorffi",
                                                        zmays = "Z. mays",
                                                        stuber = "S. tuberosum")) %>%
  mutate(species = dplyr::recode(species, atricho = "A. trichopoda",
                                 osat = "O. sativa",
                                 chlamy = "C. reinhardtii",
                                 pppatens = "P. patens",
                                 esals = "E. salsugineum",
                                 slyco = "S. lycopersicum", 
                                 athal = "A. thaliana", 
                                 smoel = "S. moellendorffi",
                                 zmays = "Z. mays",
                                 stuber = "S. tuberosum")) %>%
  ggplot(aes(y=value, x=species)) + 
  geom_bar(aes(fill=variable), color = "black", stat="identity") +
  theme_classic() + #scale_fill_manual(values=c("#C9E4CA", "#ed7a53")) +
 scale_fill_manual(values=c("#FDEF8D","#769FAE")) + 
  geom_text(data = perc_old, aes(label = correct,y = Total , x = species, vjust=-0.5), size=7) +
  labs(x="", y="# of predicted lncRNAs") + 
  guides(fill=guide_legend(title="")) + theme(legend.position = c(0.2, 0.9), 
                                              axis.text.x = element_text(face = "italic", angle = 45, hjust = 1 ),text = element_text(size=30)) +
  scale_y_continuous(limits = c(0, 8000))
ggsave("~/Documents/jobapps/figeys/talk/figs/perc_bar.pdf",
       plot=perc_bar, width = 15, heigh = 8.6, units = "in")


## again, different iterations of the same plot. some for presentations, some for paper.




