library(funrar)
library(ggplot2)
library(openxlsx)
library(vegan)
library(phytools)
library(picante)
library(geiger)
library(GUniFrac)
library(patchwork)

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
  axis.title.x = element_text(colour='black', size=13),
  axis.title.y = element_text(colour='black', size=13),
  axis.text = element_text(colour='black',size=11),
  plot.background = element_blank(), 
  plot.tag = element_text(size = 14, face = "bold")) 

###########################  Greenhouse experiment #############################
# 加载温室实验数据
green_otu = read.xlsx("Greenhouse_data_asv.xlsx", sheet = "Overall_otu", colNames = T, rowNames = T)
green_otu[1:6,1:6]
green_otu = green_otu[,-c(1)]

## load group data
Green_group = read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group-3-28", colNames = T, rowNames = T)
Green_group$Sample_ID = rownames(Green_group)
colnames(Green_group)
Green_group = Green_group[colnames(green_otu), ]
rownames(t(green_otu)) %in% rownames(Green_group)

## 群落丰富度计算
green_richness <- as.data.frame(specnumber(t(green_otu)))
green_richness$Sample_ID = rownames(green_richness); colnames(green_richness)[1] = "SR"
Green_group = Green_group %>% left_join(green_richness)
rownames(Green_group) = Green_group$Sample_ID

## Richness
green_SR_mod = lm(SR ~ Origin/Species, data = Green_group)
shapiro.test(residuals(green_SR_mod))
anova(green_SR_mod)
Richness_lm_sum = as.data.frame(anova(green_SR_mod))
Richness_lm_sum$`p.adj`=p.adjust(Richness_lm_sum$`Pr(>F)`, method = "BH")


## 群落相异性计算
green_hel_no = t(green_otu)/9690
rowSums(green_hel_no)
Bray_dist_green_no <- vegdist(green_hel_no, method = 'bray')

## perMANOVA
Green_group = Green_group[rownames(as.matrix(Bray_dist_green_no)),]

set.seed(1234)
green_perMAONA2 = GUniFrac::adonis3(green_hel_no ~ Origin/Species, method = "bray", by = "margin", 
                                    data = Green_group, permutations = 9999)
green_perMAONA2$aov.tab
permanova_green = as.data.frame(green_perMAONA2$aov.tab)
permanova_green$R2 = round(permanova_green$R2,3)
permanova_green$`F` = round(permanova_green$`F`,2)
permanova_green$`q-vaules`=p.adjust(permanova_green$`Pr(>F)`, method = "BH")
print(permanova_green)

### 功能性状与真菌组成之间关系
factors_names = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF",
                  "Funct_PC1", "Funct_PC2", "Phylo_PC1", "Phylo_PC2")

Green_Bray_Curtis_mantel = NULL
#i = "Hmax"
for (i in factors_names) {
  all_otu_TEST = Green_group[,c(i, "Sample_ID")]
  all_otu_TEST = all_otu_TEST[rownames(as.matrix(Bray_dist_green_no)), ]
  all_otu_TEST$Sample_ID = NULL
  #traits_dis <- vegdist(all_otu_TEST, method = 'euclidean')  
  factors_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE) 
  ##Bray-Curtis
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(as.dist(factors_dis), Bray_dist_green_no, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result = data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,"Fungal composition","Factors",'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result) = c("Mantel_R","P_value","from","to","Group","dist_type")
  Green_Bray_Curtis_mantel = rbind(Green_Bray_Curtis_mantel,mantel_Bray_Curtis_result)
}

# p.adjust
Green_Bray_Curtis_mantel$p.adj = p.adjust(Green_Bray_Curtis_mantel$P_value, method = "BH")
print(Green_Bray_Curtis_mantel)


##############################  Field experiment ###############################
# Soil sample grouping information
Field_group = read.xlsx("field_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$RS = sqrt(Field_group$RS)
Field_group$SRL = log10(Field_group$SRL)
shapiro.test(sqrt(Field_group$Wcont))
shapiro.test(sqrt(Field_group$Soil_N))
shapiro.test(sqrt(Field_group$Soil_ph))
Field_group$Wcont = sqrt(Field_group$Wcont)
Field_group$Soil_N = sqrt(Field_group$Soil_N)
Field_group$Years = as.factor(Field_group$Years)
Field_group$Site = as.factor(Field_group$Site)
unique(Field_group$Species)
Field_group$Site = factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin = factor(Field_group$Origin, levels = c("Native","Exotic"))

# notes: I have completed the above work, so I directly load the completed file
fungi_Flattening = read.xlsx("fungi_Flattening.xlsx", sheet = "field_flattening", rowNames = T, colNames = T)
fungi_Flattening = fungi_Flattening[,-c(1)]
fungi_Flattening[1:6, 1:6]
colSums(fungi_Flattening)

# Consider normalizing your data
Field_group = Field_group[colnames(fungi_Flattening), ] 
colnames(Field_group)
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")]))

Field_group_scale = Field_group
colnames(Field_group_scale)
Field_group_scale[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")] = 
  scale(Field_group_scale[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")])


############################## Fungal richness #################################
field_richness <- as.data.frame(specnumber(t(fungi_Flattening)))
colnames(field_richness) = "SR"
field_richness$Sample_ID = rownames(field_richness)
Field_group = Field_group %>% left_join(field_richness)
rownames(Field_group) = Field_group$Sample_ID

Fungal_SR_mod = lm(SR ~ Years + Site + Origin/Species + Years:Site + 
                     Origin:Years + Origin:Site + Origin:Years:Site + 
                     (Origin/Species):Years + (Origin/Species):Site, data = Field_group)
anova(Fungal_SR_mod)
SR_mod_anova <- as.data.frame(anova(Fungal_SR_mod))
SR_mod_anova$p_adj <- round(p.adjust(SR_mod_anova$`Pr(>F)`, method = "BH"), 3) 
SR_mod_anova[, 2:4] <- round(SR_mod_anova[, 2:4], digits = 2)
SR_mod_anova$`Pr(>F)` <- round(SR_mod_anova$`Pr(>F)`, 3)
print(SR_mod_anova) 


############################# Fungal composition ###############################
field_hel_no = t(fungi_Flattening)/9477
rowSums(field_hel_no)
Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')

set.seed(1234)
permanova_mod = GUniFrac::adonis3(field_hel_no ~ Years + Site + Origin/Species + Years:Site + 
                                    Origin:Years + Origin:Site + Origin:Years:Site + 
                                    (Origin/Species):Years + (Origin/Species):Site, method = "bray", by = "margin",
                                  data = Field_group_scale, permutations = 9999)

permanova_mod_sum = as.data.frame(permanova_mod$aov.tab)
permanova_mod_sum$R2 = round(permanova_mod_sum$R2,3)
permanova_mod_sum$df = paste0(permanova_mod_sum$Df, ",", permanova_mod_sum$Df[length(permanova_mod_sum$Df)-1])
permanova_mod_sum$`F` = round(permanova_mod_sum$F.Model,2)
permanova_mod_sum$`p.adj`= round(p.adjust(permanova_mod_sum$`Pr(>F)`, method = "BH"), 3)
print(permanova_mod_sum)

################################################################################
#### 单个性状与真菌组成之间的关系mantel_test
colnames(Field_group)
factors_names = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF",
                  "Funct_PC1", "Funct_PC2", "Phylo_PC1", "Phylo_PC2")

Field_Bray_Curtis_mantel = NULL
#i = "Hmax"
for (i in factors_names) {
  all_otu_TEST = Field_group[,c(i, "Sample_ID")]
  all_otu_TEST = all_otu_TEST[rownames(as.matrix(Bray_dist_field_no)), ]
  all_otu_TEST$Sample_ID = NULL
  #traits_dis <- vegdist(all_otu_TEST, method = 'euclidean')  
  factors_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE) 
  ##Bray-Curtis
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(as.dist(factors_dis), Bray_dist_field_no, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result = data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,"Fungal composition","Factors",'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result) = c("Mantel_R","P_value","from","to","Group","dist_type")
  Field_Bray_Curtis_mantel = rbind(Field_Bray_Curtis_mantel,mantel_Bray_Curtis_result)
}

# p.adjust
Field_Bray_Curtis_mantel$p.adj = p.adjust(Field_Bray_Curtis_mantel$P_value, method = "BH")
print(Field_Bray_Curtis_mantel)
