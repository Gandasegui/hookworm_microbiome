# Alpha diversity

```R
#Plotting refractary ASVs
otu_tab <- t(abundances(pseq))
# ASVs <- vegan::rarecurve(otu_tab, 
#                          step = 50, label = FALSE, 
#                          sample = min(rowSums(otu_tab), 
#                                       col = "blue", cex = 0.6))

plot_taxa_prevalence(pseq, "Phylum")

#rarefying at 1900 keeps all the samples, so go for it
pseq.rar <- rarefy_even_depth(pseq, rngseed=1, sample.size=1900, replace=F)

#Let's extract relevant index for alpha diversity:
#Diversity - Shannon
#Richness - Chao1
#Dominance/Evenness - simpson
#Rarity and low abundance

div_tab <- microbiome::alpha(pseq.rar, index = "all") %>%
  dplyr::select('diversity_shannon', 'chao1', 'diversity_gini_simpson', 
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
div_analysis$cured_all <-  str_replace(div_analysis$cured_all, '1', 'cured')
div_analysis$cured_all <-  str_replace(div_analysis$cured_all, '0', 'not_cured')

#LEt's do some plots
#Shannon
box_sha_1 <- div_analysis %>%
  select(diversity_shannon, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = diversity_shannon, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Shannon') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Chao1
box_chao_1 <- div_analysis %>%
  select(chao1, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = chao1, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Chao1') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Simpson
box_sim_1 <- div_analysis %>%
  select(diversity_gini_simpson, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = diversity_gini_simpson, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Simpson') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#PD
box_PD_1 <- div_analysis %>%
  select(PD, cured_all) %>%
  na.omit() %>%
  ggplot(aes(x = cured_all, y = PD, fill = cured_all)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('PD') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Let's arrange
alpha_div <- ggarrange(box_sha_1, box_chao_1, box_sim_1, box_PD_1, ncol = 4,
          labels = c('a', 'b', 'c', 'd'))
ggsave('alpha_div_cured.tiff', height = 4, width = 8)

#No diffewrences in diversity
```
# Beta diversity

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
beta_div_unw$cured_all <-  str_replace(beta_div_unw$cured_all, '1', 'cured')
beta_div_unw$cured_all <-  str_replace(beta_div_unw$cured_all, '0', 'not cured')

poca_unw_cured <- beta_div_unw %>%
  select(Axis.1, Axis.2, cured_all)%>%
  na.omit() %>%
  ggplot(., aes(x = Axis.1, y = Axis.2, col = cured_all)) + 
  geom_point(size = 3) +
  theme_light() + 
  xlab('PC1 (4.3%)') +
  ylab('PC2 (3.4%)') +
  guides(col=guide_legend(title="")) +
  scale_color_viridis_d()

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
beta_div_wt$cured_all <-  str_replace(beta_div_wt$cured_all, '1', 'cured')
beta_div_wt$cured_all <-  str_replace(beta_div_wt$cured_all, '0', 'not cured')

poca_wt_cured <- beta_div_wt %>%
  select(Axis.1, Axis.2, cured_all)%>%
  na.omit() %>%
  ggplot(aes(x = Axis.1, y = Axis.2, col = cured_all)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (12.3%)') +
  ylab('PC2 (10.2%)') + 
  scale_color_viridis_d()

beta_div <- ggarrange(poca_unw_cured, poca_wt_cured, labels = c('e', 'f'),
          common.legend = T, legend = "top")
```
Generating Figure1
```R
ggarrange(alpha_div, beta_div, nrow = 2, common.legend = F)
```
