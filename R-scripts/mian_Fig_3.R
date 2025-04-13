library(ggplot2)
mytheme = theme(
  legend.position = "none",
  panel.grid=element_blank(), 
  strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid"),
  strip.text.x = element_text(size = 9, color = "black"), # face = "bold.italic"
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=12),
  axis.title.y = element_text(colour='black', size=12),
  axis.text = element_text(colour='black',size=11),
  plot.background = element_blank(), 
  plot.tag = element_text(size = 13, face = "bold")) 

library(openxlsx)
library(picante)
traits_mean = read.xlsx("Greenhouse_data_group.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
colnames(traits_mean)
## 数据转化
shapiro.test(log10(traits_mean$SRL))
shapiro.test(sqrt(traits_mean$RS))
shapiro.test(log10(traits_mean$LCC))
shapiro.test(log10(traits_mean$LNC))
traits_mean$SRL = log10(traits_mean$SRL)
traits_mean$RS = sqrt(traits_mean$RS)


## 系统发育矩阵
plant_tree = read.tree("IQ_tree_plant_2025.NEWICK")
plot(plant_tree)
##删除某个枝长
#library(treeio)
to_drop<-c("Amborella_trichopoda","")
plant_tree <- drop.tip(as.phylo(plant_tree), to_drop) 

traits_mean = traits_mean[plant_tree$tip.label, ]

##
################################# 主成分分析
# 功能性状
select_traits2 = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF") 
pca <- prcomp(as.matrix(traits_mean[,select_traits2]), center = T, scale = T)
df_pca <- as.data.frame(pca$x)[,c(1,2)] 
df_pca$Species = rownames(df_pca)

PCOA_fun = df_pca %>% left_join(traits_mean[,c("Species", "Origin")], by = "Species") 
PCOA_fun$Origin = factor(PCOA_fun$Origin, levels = c("Native","Exotic"))

pca_loadings <- data.frame(Variables = rownames(pca$rotation),pca$rotation)
pca_loadings$label = c("Leaf chlorophyll", "Specific leaf area", "Leaf dry matter content", 
                       "Specific root length", "Fine-to-total root mass", "Root mass fraction")

PCOA_fun$Origin = factor(PCOA_fun$Origin, levels = c("Native","Exotic"))

library(ggrepel)
ggplot(data=PCOA_fun, aes(x = PC1, y = PC2))+
  geom_point(aes(shape = Origin, color = Origin, fill = Origin), size=2.2, show.legend = T)+
  geom_segment(data = pca_loadings, mapping = aes(x = 0, y = 0, xend = (PC1*4), yend = (PC2*4)), 
               color = "black", alpha=0.8) +
  annotate("text", x = (pca_loadings$PC1*4.5), y = (pca_loadings$PC2*4.5),label = pca_loadings$Variables, size = 3.5) +
  labs(x = paste("PCA1 (", sprintf("%.1f", summary(pca)$importance[2, 1] * 100), "%)", sep = ""),
       y = paste("PCA2 (", sprintf("%.1f", summary(pca)$importance[2, 2] * 100), "%)", sep = ""), tag = "b") +
  geom_hline(yintercept=0,linetype=2)+ 
  geom_vline(xintercept=0,linetype=2)+
  theme_classic() + mytheme + 
  theme(legend.position=c(0.85,0.95)) + 
  scale_shape_manual(values = c(16,15)) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) -> Fig_2b; Fig_2b

library(ggExtra)
Fig_2bb <- ggMarginal(Fig_2b, type = "boxplot", groupColour = F, 
                      groupFill = TRUE, alpha = 1, color = "black", size = 8) ; Fig_2bb

################################################################################
## 系统发育矩阵
plant_tree = read.tree("IQ_tree_plant_2025.NEWICK")
##删除某个枝长
#library(treeio)
to_drop<-c("Amborella_trichopoda","")
plant_tree <- drop.tip(as.phylo(plant_tree), to_drop) 

plant_dist <- cophenetic.phylo(plant_tree)
#plant_dist <- cophenetic(plant_tree)

pca <- prcomp(plant_dist)

df_pca <- as.data.frame(pca$x)[,c(1,2)] 
df_pca$Species = rownames(df_pca)

#write.csv(df_pca,"phylo_df_pca.csv")

PCOA_phylo = df_pca %>% left_join(traits_mean[,c("Species", "Origin")], by = "Species") 
PCOA_phylo$Origin = factor(PCOA_phylo$Origin, levels = c("Native","Exotic"))

library(ggrepel)
ggplot(data=PCOA_phylo, aes(x = PC1, y = PC2))+
  geom_point(aes(shape = Origin, color = Origin, fill = Origin), size=2.2, show.legend = T)+
  labs(x = paste("PCA1 (", sprintf("%.1f", summary(pca)$importance[2, 1] * 100), "%)", sep = ""),
       y = paste("PCA2 (", sprintf("%.1f", summary(pca)$importance[2, 2] * 100), "%)", sep = ""), tag = "c") +
  geom_hline(yintercept=0,linetype=2)+ 
  geom_vline(xintercept=0,linetype=2)+
  theme_classic() + mytheme + 
  theme(legend.position=c(0.12,0.85)) + 
  scale_shape_manual(values = c(16,15)) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) -> Fig_2c; Fig_2c

library(ggExtra)
Fig_2cc <- ggMarginal(Fig_2c, type = "boxplot", groupColour = F, 
                      groupFill = TRUE, alpha = 1, color = "black", size = 8) ; Fig_2cc


### 性状差异比较pgls
library(caper)
library(dplyr)
PCOA_phylo2 = PCOA_phylo
PCOA_fun2 = PCOA_fun
colnames(PCOA_fun2)[1:2] = c("Funct_PC1", "Funct_PC2")
colnames(PCOA_phylo2)[1:2] = c("Phylo_PC1", "Phylo_PC2")
traits_mean_pca = traits_mean %>% left_join(PCOA_fun2) %>% left_join(PCOA_phylo2)

vector <- c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF","Funct_PC1", "Funct_PC2","Phylo_PC1", "Phylo_PC2")

pgls_final <- NULL
for (i in vector) {
  traits_all1 <- traits_mean_pca[,c(i, "Origin", "Species")]
  pgls_data <- comparative.data(phy = plant_tree, data = traits_all1, 
                                names.col = Species, vcv = TRUE, 
                                na.omit = FALSE, warn.dropped = TRUE)
  colnames(pgls_data$data)[1] <- "Traits"
  model.pgls <- pgls(scale(Traits) ~ Origin, data = pgls_data, lambda = "ML")
  pgls_coef <- summary(model.pgls)$coef[2]
  pgls_se <- summary(model.pgls)$coef[2,2]
  pgls_pvalue <- summary(model.pgls)$coef[2,4]
  pgls_summry <- data.frame(Traits = i, coef = pgls_coef, se = pgls_se, p = pgls_pvalue)
  pgls_final <- rbind(pgls_final, pgls_summry)
}

i = "Funct_PC1"
anova(model.pgls)
########################## 性状比较可视化 ######################################
traits_mean_long <- traits_mean_pca %>%
  tidyr::pivot_longer(cols = c(Chol, SLA, LDMC, SRL, FRR, RMF),  # 选择需要转换的列
                      names_to = "Traits",       # 新列的名称
                      values_to = "Traits_values")  # 新值的列名


library(gghalves)
traits_mean_long$Origin = factor(traits_mean_long$Origin, levels = c("Native", "Exotic"))
traits_mean_long$Traits = factor(traits_mean_long$Traits, levels = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF"))

## 更换标签
Traits_label = c(Chol = "Leaf chlorophyll (SPAD)", SLA = "Specific leaf area (cm2 g-1)", LDMC = "Leaf dry matter content (g g-1)",
                 SRL = "Specific root length (cm g-1, log10)", FRR = "Fine-to-total root mass (g g-1)", RMF = "Root mass fraction (g g-1)")

ggplot(traits_mean_long, aes(x = Origin,y = Traits_values, fill = Origin, shape = Origin))+
  geom_half_violin(position=position_nudge(x=0.15, y=0),side='R',adjust=1.2,trim=T,color=NA,alpha=0.8) +
  geom_boxplot(width = 0.2, alpha = 1, outliers = FALSE) + 
  geom_point(aes(x = Origin, y = Traits_values, fill = Origin), color = "black",  
             position = position_jitter(width = 0.08), size = 1.5,alpha = 0.6) + 
  labs(x = '', y = 'Traits value', tag = "a") + 
  theme_bw() + mytheme + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 35)) + 
  scale_shape_manual(values = c(21, 22)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) + 
  facet_wrap(Traits ~., ncol = 2, nrow = 5, scales = "free_y",
             labeller = labeller(Traits = Traits_label)) -> Fig_2a; Fig_2a







library(patchwork)
(Fig_2a|(Fig_2b/Fig_2c)) + plot_layout(widths = c(0.5, 0.5))

################################################################################






