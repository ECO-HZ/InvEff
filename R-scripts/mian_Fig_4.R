library(ggplot2)
library(ggtext) # version 0.1.2
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
  plot.title = element_textbox(
    size = 14, color = "black", fill = "grey90",
    box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
  plot.tag = element_text(size = 14, face = "bold")) 

library(openxlsx)
library(adespatial)
library(vegan)
library(phytools)
library(picante)
library(geiger)
# Soil sample grouping information
Field_group = read.xlsx("field_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$RS = sqrt(Field_group$RS)
Field_group$SRL = log10(Field_group$SRL)
Field_group$Wcont = sqrt(Field_group$Wcont)
Field_group$Soil_N = sqrt(Field_group$Soil_N)
Field_group$Years = as.factor(Field_group$Years)
Field_group$Site = as.factor(Field_group$Site)
Field_group = Field_group[colnames(fungi_Flattening), ] 
unique(Field_group$Species)
Field_group$Site = factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin = factor(Field_group$Origin, levels = c("Native","Exotic"))

# notes: I have completed the above work, so I directly load the completed file
fungi_Flattening = read.xlsx("fungi_Flattening.xlsx", sheet = "field_flattening", rowNames = T, colNames = T)
fungi_Flattening = fungi_Flattening[,-c(1)]
fungi_Flattening[1:6, 1:6]
colSums(fungi_Flattening)

# Consider normalizing your data
colnames(Field_group)
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")]))

Field_group_scale = Field_group
colnames(Field_group_scale)
Field_group_scale[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")] = 
  scale(Field_group_scale[c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec","Hmax","Chol","LA","SLA","LDMC","SRL","FRR","RS","RMF","AGB","BGB","TB","LCC","LNC","LCN","Fun_PC1","Fun_PC2","Phylo_PC1","Phylo_PC2","Fun_RC1","Fun_RC2","Funct_PC1","Funct_PC2","Funct_RC1","Funct_RC2")])

## not hellinger
field_hel_no = t(fungi_Flattening)/9477
rowSums(field_hel_no)
Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')

################################################################################
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

## 群落丰度度计算
green_richness <- as.data.frame(specnumber(t(green_otu)))
green_richness$Sample_ID = rownames(green_richness); colnames(green_richness)[1] = "SR"
Green_group = Green_group %>% left_join(green_richness)
rownames(Green_group) = Green_group$Sample_ID

## 群落相异性计算
green_hel_no = t(green_otu)/9690
rowSums(green_hel_no)
Bray_dist_green_no <- vegdist(green_hel_no, method = 'bray')


################################################################################
## 计算物种bray距离均值
Green_dist_data <- reshape2::melt(as.matrix(Bray_dist_green_no), varnames = c("Sample_ID_A", "Sample_ID_B"),
                                  value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] = "Sample_ID"
### 
Green_dist_data = Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] = c("Sample_ID2","Sample_ID")
Green_dist_data = Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
Green_dist_re = Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
# 将长数据还原为对称矩阵
Bray_dist_green_mean <- reshape2::dcast(Green_dist_re, Species.x  ~ Species.y , value.var = "dist")
rownames(Bray_dist_green_mean) = Bray_dist_green_mean$Species.x
Bray_dist_green_mean = Bray_dist_green_mean[,-1]
diag(Bray_dist_green_mean) = 0

################################################################################
## deming 回归 ### 采用Bray距离进行评估

standr = function(x){(x-min(x))/(max(x)-min(x))} 

## pisewise - Native & Native
Years = unique(Field_group_scale$Years)
Site = unique(Field_group_scale$Site)
β_Bray_Native = NULL

for(i in Years){
  for(ii in Site) {
    ## Field
    select_group = subset(Field_group, Years == i & Site == ii)
    ##
    native_sample = subset(select_group, Origin == "Native")$Sample_ID
    native_latin = subset(select_group, Origin == "Native")$Species
    ## Field
    select_dist_field = as.matrix(Bray_dist_field_no)[native_sample, native_sample]
    select_dist_field = standr(select_dist_field)
    Field_mean = mean(select_dist_field)
    Field_sd = sd(select_dist_field)
    Field_se = Field_sd / sqrt(length(select_dist_field))  # 计算标准误差 (SE)
    ## Greenhouse 
    select_dist_green = as.matrix(Bray_dist_green_mean)[native_latin, native_latin]
    select_dist_green <- standr(select_dist_green)
    Green_mean = mean(select_dist_green)
    Green_sd = sd(select_dist_green)
    Green_se = Green_sd / sqrt(length(select_dist_green))  # 计算标准误差 (SE)
    
    ## 合并数据
    dist_data_all = data.frame(Years = i , Site = ii, 
                               Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                               Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se)
    β_Bray_Native = rbind(β_Bray_Native, dist_data_all)
  }
}
# View(β_Bray_Native)

library(deming)
fit <- deming::deming(Field_mean ~ Green_mean, ystd=Field_se, xstd=Green_se, data=β_Bray_Native)
print(fit)
β_Bray_Native$Site = factor(β_Bray_Native$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept=0,slope=1, color = "#8B0000", linetype = 1, size=1)+
  geom_point(β_Bray_Native, mapping = aes(x=Green_mean,y=Field_mean, shape = Years, fill = Site), color = "black", size=2.5)+
  geom_errorbar(data = β_Bray_Native,mapping = aes(x = Green_mean,ymax = Field_mean+Field_se, ymin=Field_mean-Field_se, color = Site),width=0,size=0.5,alpha = 1)+#
  geom_errorbarh(data = β_Bray_Native,mapping = aes(y = Field_mean,xmax=Green_mean+Green_se,xmin=Green_mean-Green_se, color = Site),height=0,size=0.5,alpha = 1)+#
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  geom_abline(intercept=-0.5805579, slope=1.9401005, linetype = 1, size=1)+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.60,0.87)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.60,0.87)) + 
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = NULL, 
       y = "Mean pair-wise Bray–Curtis dissimilarities\nestimated in the field",
       tag = "a", title = c("Between native plants")) -> mian_Fig_4a; mian_Fig_4a

## pisewise - Exotic & Exotic
Years = unique(Field_group_scale$Years)
Site = unique(Field_group_scale$Site)
β_Bray_Exotic = NULL

for(i in Years){
  for(ii in Site) {
    ## Field
    select_group = subset(Field_group, Years == i & Site == ii)
    ##
    exotic_sample = subset(select_group, Origin == "Exotic")$Sample_ID
    exotic_latin = subset(select_group, Origin == "Exotic")$Species
    ## Field
    select_dist_field = as.matrix(Bray_dist_field_no)[exotic_sample, exotic_sample]
    select_dist_field = standr(select_dist_field)
    Field_mean = mean(select_dist_field)
    Field_sd = sd(select_dist_field)
    Field_se = Field_sd / sqrt(length(select_dist_field))  # 计算标准误差 (SE)
    ## Greenhouse 
    select_dist_green = as.matrix(Bray_dist_green_mean)[exotic_latin, exotic_latin]
    select_dist_green <- standr(select_dist_green)
    Green_mean = mean(select_dist_green)
    Green_sd = sd(select_dist_green)
    Green_se = Green_sd / sqrt(length(select_dist_green))  # 计算标准误差 (SE)
    
    ## 合并数据
    dist_data_all = data.frame(Years = i , Site = ii, 
                               Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                               Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se)
    β_Bray_Exotic = rbind(β_Bray_Exotic, dist_data_all)
  }
}
# View(β_Bray_Exotic)

library(deming)
fit <- deming::deming(Field_mean ~ Green_mean, ystd=Field_se, xstd=Green_se, data=β_Bray_Exotic)
print(fit)
β_Bray_Exotic$Site = factor(β_Bray_Exotic$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept=0,slope=1, color = "#96393D", linetype = 1, size=1)+
  geom_point(β_Bray_Exotic, mapping = aes(x=Green_mean,y=Field_mean, shape = Years, fill = Site), color = "black", size=2.5)+
  geom_errorbar(data = β_Bray_Exotic,mapping = aes(x = Green_mean,ymax = Field_mean+Field_se, ymin=Field_mean-Field_se, color = Site),width=0,size=0.5,alpha = 1)+#
  geom_errorbarh(data = β_Bray_Exotic,mapping = aes(y = Field_mean,xmax=Green_mean+Green_se,xmin=Green_mean-Green_se, color = Site),height=0,size=0.5,alpha = 1)+#
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  geom_abline(intercept=5.774493, slope=-6.902086, linetype = 2, size=1)+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.60,0.87)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.60,0.87)) + 
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = "Mean pair-wise Bray–Curtis dissimilarities\nestimated in the greenhouse experiment", 
       y = NULL, tag = "b", title = c("Between exotic plants")) -> mian_Fig_4b; mian_Fig_4b


######################### pisewise - Native & Exotic ###########################
Years = unique(Field_group_scale$Years)
Site = unique(Field_group_scale$Site)
β_Bray_data = NULL

for(i in Years){
  for(ii in Site) {
    ## Field
    select_group = subset(Field_group, Years == i & Site == ii)
    ##
    native_sample = subset(select_group, Origin == "Native")$Sample_ID
    native_latin = subset(select_group, Origin == "Native")$Species
    exotic_sample = subset(select_group, Origin == "Exotic")$Sample_ID
    exotic_latin = subset(select_group, Origin == "Exotic")$Species
    ## Field
    select_dist_field = as.matrix(Bray_dist_field_no)[native_sample, exotic_sample]
    select_dist_field <- standr(select_dist_field)
    Field_mean = mean(select_dist_field)
    Field_sd = sd(select_dist_field)
    Field_se = Field_sd / sqrt(length(select_dist_field))  # 计算标准误差 (SE)
    ## Greenhouse 
    select_dist_green = as.matrix(Bray_dist_green_mean)[native_latin, exotic_latin]
    select_dist_green <- standr(select_dist_green)
    Green_mean = mean(select_dist_green)
    Green_sd = sd(select_dist_green)
    Green_se = Green_sd / sqrt(length(select_dist_green))  # 计算标准误差 (SE)
    
    ## 合并数据
    dist_data_all = data.frame(Years = i , Site = ii, 
                               Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                               Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se)
    β_Bray_data = rbind(β_Bray_data, dist_data_all)
  }
}


library(deming)
fit <- deming::deming(Field_mean ~ Green_mean, ystd=Field_se, xstd=Green_se, data=β_Bray_data)
print(fit)
β_Bray_data$Site = factor(β_Bray_data$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept=0,slope=1, color = "#96393D", linetype = 1, size=1)+
  geom_point(β_Bray_data, mapping = aes(x=Green_mean,y=Field_mean, shape = Years, fill = Site), color = "black", size=2.5)+
  geom_errorbar(data = β_Bray_data,mapping = aes(x = Green_mean,ymax = Field_mean+Field_se, ymin=Field_mean-Field_se, color = Site),width=0,size=0.5,alpha = 1)+#
  geom_errorbarh(data = β_Bray_data,mapping = aes(y = Field_mean,xmax=Green_mean+Green_se,xmin=Green_mean-Green_se, color = Site),height=0,size=0.5,alpha = 1)+#
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  geom_abline(intercept=2.212175, slope=-3.098372, linetype = 2, size=1)+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.25,0.88)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.25,0.88)) + 
  theme_bw() + mytheme + theme(legend.position = c(0.85,0.45)) + 
  labs(x = NULL, 
       y = NULL, tag = "c", title = c("Between native and exotic plants")) -> mian_Fig_4c; mian_Fig_4c



mian_Fig_4a|mian_Fig_4b|mian_Fig_4c

