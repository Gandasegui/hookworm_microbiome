# Rscript to peform microbiome analysis

```R
library(tidyverse)
library(microbiome)
library(phyloseq)
library(knitr)
library(metagMisc)
library(picante)
library(data.table)
library(ggpubr)
library(vegan)
library(umap)
```

### Alpha diversity

```R
#Plotting refractary ASVs
otu_tab <- t(abundances(pseq))
ASVs <- vegan::rarecurve(otu_tab, 
                         step = 50, label = FALSE, 
                         sample = min(rowSums(otu_tab), 
                                      col = "blue", cex = 0.6))

p.rar <- plot_taxa_prevalence(pseq, "Phylum")

#rarefying at 1100 keeps all the samples, so go for it
pseq.rar <- rarefy_even_depth(pseq, rngseed=1, sample.size=1100, replace=F)

#Let's extract relevant index for alpha diversity:
#Diversity - Shannon
#Richness - Chao1
#Dominance/Evenness - simpson
#Rarity and low abundance

div_tab <- microbiome::alpha(pseq.rar, index = "all") %>%
  select('diversity_shannon', 'chao1', 'diversity_gini_simpson', 
         'diversity_inverse_simpson', 'observed') %>%
  rownames_to_column(var = 'sample-id')

#Now, let's calculate faith-PD
pseq.asvtab <- as.data.frame(pseq.rar@otu_table)
pseq.tree <- pseq.rar@phy_tree #check it is rooted

div_pd <- pd(t(pseq.asvtab), pseq.tree, include.root=T) %>%
  rownames_to_column(var = 'sample-id')
#t(ou_table) transposes the table for use in picante
# and the tree file comes from the first code chunck we used to read tree file
#(see making a phyloseq object section).
print(div_pd)

#Now, we modify the metadata for joining
div_meta <- rownames_to_column(metadata, var = 'sample-id')

#And joining
div_analysis <- left_join(div_meta, div_tab, by = 'sample-id') %>%
  left_join(., div_pd, by = 'sample-id')

div_analysis$cured_all <- as.factor(div_analysis$cured_all)
div_analysis$cured_all <-  str_replace(div_analysis$cured_all, '0', 'Not Cured')
div_analysis$cured_all <-  str_replace(div_analysis$cured_all, '1', 'Cured')

#LEt's do some plots
#Shannon
box_sha_1 <- div_analysis %>%
  select(diversity_shannon, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = diversity_shannon, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Shannon') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none") +
  labs(x ='')

#Chao1
box_chao_1 <- div_analysis %>%
  select(chao1, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = chao1, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Chao1') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")  +
  labs(x ='')

#Simpson
box_sim_1 <- div_analysis %>%
  select(diversity_gini_simpson, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = diversity_gini_simpson, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Simpson') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")  +
  labs(x ='')

#PD
box_PD_1 <- div_analysis %>%
  select(PD, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = PD, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('PD') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")  +
  labs(x ='')

#Let's arrange
alpha_div <- ggarrange(box_sha_1, box_chao_1, box_sim_1, box_PD_1, ncol = 4,
          labels = c('a'))
#No differences in diversity
```
### Beta diversity

```R
#Unweighted Unifrac
ordu.unwt.uni <- ordinate(pseq.rar, "PCoA", "unifrac", weighted=F)
barplot(ordu.unwt.uni$values$Eigenvalues[1:10])
barplot(ordu.unwt.uni$values$Relative_eig[1:10])

#Plotting to know the variance explained
unwt.unifrac <- plot_ordination(pseq.rar, ordu.unwt.uni) +
  ggtitle("Unweighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer()


eigen_vec_unw <- rownames_to_column(data.frame(ordu.unwt.uni$vectors), var = 'sample-id')
beta_div_unw <- left_join(div_meta, eigen_vec_unw, by = 'sample-id')
beta_div_unw$cured_all <- as.factor(beta_div_unw$cured_all)
beta_div_unw$cured_all <-  str_replace(beta_div_unw$cured_all, '1', 'Cured')
beta_div_unw$cured_all <-  str_replace(beta_div_unw$cured_all, '0', 'Not cured')

poca_unw_cured <- beta_div_unw %>%
  select(Axis.1, Axis.2, cured_all)%>%
  na.omit() %>%
  ggplot(., aes(x = Axis.1, y = Axis.2, col = cured_all)) + 
  geom_point(size = 3) +
  theme_light() + 
  xlab('PC1 (8.2%)') +
  ylab('PC2 (7.8%)') +
  guides(col=guide_legend(title=""))

#Weighted Unifrac
ordu.wt.uni <- ordinate(pseq.rar, "PCoA", "unifrac", weighted=T)
barplot(ordu.wt.uni$values$Eigenvalues[1:10])
barplot(ordu.wt.uni$values$Relative_eig[1:10])

#Plotting to know the variance explained
wt.unifrac <- plot_ordination(pseq.rar, ordu.wt.uni) +
  ggtitle("Weighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer()


eigen_vec_wt <- rownames_to_column(data.frame(ordu.wt.uni$vectors), var = 'sample-id')
beta_div_wt <- left_join(div_meta, eigen_vec_wt, by = 'sample-id')
beta_div_wt$cured_all <- as.factor(beta_div_wt$cured_all)
beta_div_wt$cured_all <-  str_replace(beta_div_wt$cured_all, '1', 'Cured')
beta_div_wt$cured_all <-  str_replace(beta_div_wt$cured_all, '0', 'Not cured')

poca_wt_cured <- beta_div_wt %>%
  select(Axis.1, Axis.2, cured_all)%>%
  na.omit() %>%
  ggplot(aes(x = Axis.1, y = Axis.2, col = cured_all)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (49.4%)') +
  ylab('PC2 (18%)') 

beta_div <- ggarrange(poca_unw_cured, poca_wt_cured, labels = c('b'),
          common.legend = T, legend = 'bottom')
ggarrange(alpha_div, beta_div, nrow = 2, heights = c(1, 1.2))
ggsave('Figures/Fig1_div.jpg', width = 8, height = 7)
ggsave('Figures/Fig1_div.tiff', width = 8, height = 7)
ggsave('Figures/Fig1_div.pdf', width = 8, height = 7)
```

### Compositional analysis - KW test

```R
#rename the variable
pseq@sam_data[["cured_all"]] <- str_replace(pseq@sam_data[["cured_all"]], '1', 'cured')
pseq@sam_data[["cured_all"]] <- str_replace(pseq@sam_data[["cured_all"]], '0', 'not_cured')

# Make sure we use functions from correct package
transform <- microbiome::transform

#PHYLUM

# Merge rare taxa to speed up examples
pseq.comp <- transform(pseq, "compositional")
pseq.comp.phy <- aggregate_rare(pseq.comp, level = "Phylum", detection = 1/100, prevalence = 10/100)

#Now, lets plot the composition sorted by LCN-2
plot_composition(pseq.comp.phy) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic"))

#Now, let's generate a new tibble for statistical analysis and representation
comp_phy_df <- t(as.data.frame(otu_table(pseq.comp.phy))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

comp_phy_df <- rownames_to_column(metadata, var = 'sample-id') %>%
  select('sample-id', 'cured_all') %>%
  left_join(comp_phy_df, ., by = 'sample-id')

#Kw
pseq_all_phy <- colnames(comp_phy_df)[-1]
pseq_all_phy <- pseq_all_phy[-8]


comp_kw_phy = data_frame(phylum = character(), pval = double())
for (i in pseq_all_phy){
  # Creating new variables makes the code clearer
  x <- comp_phy_df %>% pull(i)
  kw <- kruskal.test(x ~ comp_phy_df$cured_all)
  pval <- kw$p.value
  comp_kw_phy <- comp_kw_phy %>% add_row('phylum' = i, 'pval' = pval)
}
#Let's adjust by bonferroni
p_adj_bon <- comp_kw_phy %>% pull(pval)
p_adj_bon <- p.adjust(p_adj_bon, method = 'BH')
comp_kw_phy$padjusted <- p_adj_bon
#nothing significant

#FAMILY

pseq.comp.fam <- aggregate_rare(pseq.comp, level = "Family", detection = 1/100, prevalence = 10/100)

plot_composition(pseq.comp.fam,
                 sample.sort = "cured_all",
                 transform = "compositional") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic"))

#Now, let's generate a new tibble for statistical analysis and representation
comp_fam_df <- t(as.data.frame(otu_table(pseq.comp.fam))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

colnames(comp_fam_df) <- gsub("[", "", colnames(comp_fam_df), fixed = TRUE)
colnames(comp_fam_df) <- gsub("]", "", colnames(comp_fam_df), fixed = TRUE)

#Now, let's generate a new tibble for statistical analysis and representation
comp_fam_df <- t(as.data.frame(otu_table(pseq.comp.fam))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

comp_fam_df <- rownames_to_column(metadata, var = 'sample-id') %>%
  select('sample-id', 'cured_all') %>%
  left_join(comp_fam_df, ., by = 'sample-id')

#Kw
pseq_all_fam <- colnames(comp_fam_df)[-1]
pseq_all_fam <- pseq_all_fam[-34]


comp_kw_fam = data_frame(family = character(), pval = double())
for (i in pseq_all_fam){
  # Creating new variables makes the code clearer
  x <- comp_fam_df %>% pull(i)
  kw <- kruskal.test(x ~ comp_fam_df$cured_all)
  pval <- kw$p.value
  comp_kw_fam <- comp_kw_fam %>% add_row('family' = i, 'pval' = pval)
}
#Let's adjust by bonferroni
p_adj_bon <- comp_kw_fam %>% pull(pval)
p_adj_bon <- p.adjust(p_adj_bon, method = 'bonferroni')
comp_kw_fam$padjusted <- p_adj_bon
#nothing significant

#GENUS

pseq.comp.gen <- aggregate_rare(pseq.comp, level = "Genus", detection = 1/100, prevalence = 10/100)

plot_composition(pseq.comp.gen,
                 sample.sort = "cured_all",
                 transform = "compositional") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic"))

#Now, let's generate a new tibble for statistical analysis and representation
comp_gen_df <- t(as.data.frame(otu_table(pseq.comp.gen))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

colnames(comp_gen_df) <- gsub("[", "", colnames(comp_gen_df), fixed = TRUE)
colnames(comp_gen_df) <- gsub("]", "", colnames(comp_gen_df), fixed = TRUE)

comp_gen_df <- rownames_to_column(metadata, var = 'sample-id') %>%
  select('sample-id', 'cured_all') %>%
  left_join(comp_gen_df, ., by = 'sample-id')

#Kw
pseq_all_gen <- colnames(comp_gen_df)[-1]
pseq_all_gen <- pseq_all_gen[-58]


comp_kw_gen = data_frame(genus = character(), pval = double())
for (i in pseq_all_gen){
  # Creating new variables makes the code clearer
  x <- comp_gen_df %>% pull(i)
  kw <- kruskal.test(x ~ comp_gen_df$cured_all)
  pval <- kw$p.value
  comp_kw_gen <- comp_kw_gen %>% add_row('genus' = i, 'pval' = pval)
}
#Let's adjust by bonferroni
p_adj_bon <- comp_kw_gen %>% pull(pval)
p_adj_bon <- p.adjust(p_adj_bon, method = 'bonferroni')
comp_kw_gen$padjusted <- p_adj_bon
#nothing significant
```

### Compositional analysis - ALDE-x2

```R
library(mia)
library('ALDEx2')
#Let's tranform the data
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)

#Phylum

#Agglomerate by genus and subset by prevalence
tse_aldex_phy <- subsetByPrevalentTaxa(tse,
                                       rank = "Phylum",
                                       prevalence = 10/100)
tse_aldex_phy <- mia::transformCounts(tse_aldex_phy, method = "relabundance")

#R un the test
aldex_phy <- aldex(assay(tse_aldex_phy), tse_aldex_phy$cured_all, denom = 'all',
                   test = 'kw')


#Family

#Agglomerate by genus and subset by prevalence
tse_aldex_fam <- subsetByPrevalentTaxa(tse,
                                       rank = "Family",
                                       prevalence = 10/100)
tse_aldex_fam <- mia::transformCounts(tse_aldex_fam, method = "relabundance")

#R un the test
aldex_fam <- aldex(assay(tse_aldex_fam), tse_aldex_fam$cured_all, denom = 'all',
                   test = 'kw')

# Genus

# Agglomerate by genus and subset by prevalence
tse_aldex_gen <- subsetByPrevalentTaxa(tse,
                                       rank = "Genus",
                                       prevalence = 10/100)
tse_aldex_gen <- mia::transformCounts(tse_aldex_gen, method = "relabundance")

#R un the test
aldex_gen <- aldex(assay(tse_aldex_gen), tse_aldex_gen$cured_all, denom = 'all',
                   test = 'kw') 

#Nothing significant in any level
```

### AUC plot from ML

```R
#AUC
auc_all <- read_csv('ML_AUCs/all taxon results.csv') %>%
  mutate(., level = 'All')

auc_phy <- read_csv('ML_AUCs/Phylum validatio results.csv') %>%
  mutate(., level = 'Phylum')

auc_fam <- read_csv('ML_AUCs/Family validatio results.csv') %>%
  mutate(., level = 'Family')

auc_gen <- read_csv('ML_AUCs/Genus validatio results.csv') %>%
  mutate(., level = 'Genus')

auc_combined <- rbind(auc_phy, auc_fam, auc_gen)

colnames(auc_combined) <- c('Model', 'Mean', 'sd', 'level')

#Let's do some mod for plotting

#LEt's use acro. for the models
auc_combined$Model <- str_replace(auc_combined$Model, 'BernoulliNB', 'BNB')
auc_combined$Model <- str_replace(auc_combined$Model, 'BaggingClassifier', 'BaC')
auc_combined$Model <- str_replace(auc_combined$Model, 'GaussianNB', 'GNB')
auc_combined$Model <- str_replace(auc_combined$Model, 'LogisticRegression', 'LoR')
auc_combined$Model <- str_replace(auc_combined$Model, 'QuadraticDiscriminantAnalysis', 'QDA')
auc_combined$Model <- str_replace(auc_combined$Model, 'RandomForestClassifier', 'RFC')

#Let's order the variables
auc_combined$Model <- factor(auc_combined$Model, levels = c('BNB', 'GNB', 
                                                            'BaC', 'QDA',
                                                            'LoR', 'RFC',
                                                            'XGB'))
auc_combined$level <- factor(auc_combined$level, levels = c('Phylum',  
                                                            'Family', 
                                                            'Genus'))

AUC_plot <- ggplot(auc_combined , aes(x = Model, y = Mean, color = level, fill = level)) + 
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=0.2, color = 'grey20',
                position = position_dodge(width=0.75)) +
  geom_point(size = 3, aes(color = level, fill = level),
             position = position_dodge(width=0.75)) +
  scale_color_manual(values=c('#43CD80','#E69F00', '#56B4E9', '#EEDC82'),
                     guide = guide_legend(override.aes = list(c('All', 'Phylum', 'Family', 'Genus')))) +
  geom_hline(yintercept = 0.7, colour = 'grey30', linetype = 'dashed', linewidth = .5) +
  geom_hline(yintercept = 0.5, colour = 'grey30', linetype = 'dashed', linewidth = .5) +
  ylim(0, 0.95) +
  facet_wrap(~ level, nrow = 1) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Model') +
  ylab("AUC")
```

### UMAP analysis

```R
#phylum
phy_umap_df <- comp_phy_df[, -1]
phy_umap_df <- as_tibble(scale(phy_umap_df[, -8]))
  
phy_umap <- umap(phy_umap_df)
phy_umap.val <- as_tibble(phy_umap$layout)

phy_umap.val <- comp_phy_df %>%
  dplyr::select(., 'cured_all') %>%
  cbind(phy_umap.val, .)%>%
  as_tibble()

phy_umap.val$cured_all <-  str_replace(phy_umap.val$cured_all, '1', 'Cured')
phy_umap.val$cured_all <-  str_replace(phy_umap.val$cured_all, '0', 'Not cured')

phy_umap.val <- mutate(phy_umap.val, level = 'Phylum')

#Family
fam_umap_df <- comp_fam_df[, -1]
fam_umap_df <- as_tibble(scale(fam_umap_df[, -34]))

fam_umap <- umap(fam_umap_df)
fam_umap.val <- as_tibble(fam_umap$layout)

fam_umap.val <- comp_fam_df %>%
  dplyr::select(., 'cured_all') %>%
  cbind(fam_umap.val, .)%>%
  as_tibble()

fam_umap.val$cured_all <-  str_replace(fam_umap.val$cured_all, '1', 'Cured')
fam_umap.val$cured_all <-  str_replace(fam_umap.val$cured_all, '0', 'Not cured')

fam_umap.val <- mutate(fam_umap.val, level = 'Family')

#Genus
gen_umap_df <- comp_gen_df[, -1]
gen_umap_df <- as_tibble(scale(gen_umap_df[, -58]))

gen_umap <- umap(gen_umap_df)
gen_umap.val <- as_tibble(gen_umap$layout)

gen_umap.val <- comp_gen_df %>%
  dplyr::select(., 'cured_all') %>%
  cbind(gen_umap.val, .)%>%
  as_tibble()

gen_umap.val$cured_all <-  str_replace(gen_umap.val$cured_all, '1', 'Cured')
gen_umap.val$cured_all <-  str_replace(gen_umap.val$cured_all, '0', 'Not cured')

gen_umap.val <- mutate(gen_umap.val, level = 'Genus')

#Let's megrge and plot
all_umap.val <- rbind(phy_umap.val, fam_umap.val, gen_umap.val)

all_umap.val$level <- factor(all_umap.val$level, levels = c('Phylum',  
                                                            'Family', 
                                                            'Genus'))

umap_plot <- ggplot(all_umap.val, aes(x=V1, y=V2, color = cured_all)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(color = "") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~level, nrow = 1) +
  labs(x = 'V1') +
  ylab("V2")

ggarrange(umap_plot, AUC_plot, nrow = 2, heights = c(1, 1.2), widths = c(1, 0.8),
          labels = c('a', 'b'))
ggsave('Figures/Fig2_ml.jpg', width = 8, height = 7)
ggsave('Figures/Fig2_ml.tiff', width = 8, height = 7)
ggsave('Figures/Fig2_ml.pdf', width = 8, height = 7)
```
