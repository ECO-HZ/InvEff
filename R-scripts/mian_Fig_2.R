library(ggplot2)
library(openxlsx)
library(adespatial)
library(vegan)
library(phytools)
library(picante)
library(geiger)
library(GUniFrac)
library(patchwork)
library(car)

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

################################## NMDS ########################################
set.seed(1234)
field_nmds1 <- metaMDS(Bray_dist_field_no, k = 3, trymax = 500, autotransform = TRUE, wascores = TRUE)
field_nmds1.stress <- field_nmds1$stress
field_plot_data <- data.frame(field_nmds1$point)
field_plot_data$Sample_ID <- rownames(field_plot_data)
names(field_plot_data)[1:3] <- c('NMDS1', 'NMDS2', 'NMDS3')
field_plot_data <- merge(field_plot_data, Field_group, by = 'Sample_ID', all.x = TRUE)
field_plot_data$Origin <- factor(field_plot_data$Origin, levels = c("Native","Exotic"))
field_plot_data$Site <- factor(field_plot_data$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot(field_plot_data, aes(NMDS1, NMDS2, shape = Origin, fill = Years))+
  geom_point(size = 2.2, color = "black", alpha = 1, show.legend = T) +
  scale_shape_manual(values = c(21,22)) +
  scale_color_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  scale_fill_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text=element_text(size=12))+
  theme_bw() + mytheme +
  theme(legend.position = "right") + 
  labs(x="NMDS1", y="NMDS2",tag = "a")+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> mian_Fig_2a; mian_Fig_2a

ggplot(field_plot_data, aes(NMDS1, NMDS2, shape = Origin, fill = Site))+
  geom_point(size = 2.2, color = "black",alpha = 1, show.legend = T) +
  scale_shape_manual(values = c(21,22)) +
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 12))+
  theme_bw() + mytheme +
  theme(legend.position = "right") + 
  labs(x = "NMDS1", y = "NMDS2")+
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) -> mian_Fig_2aa; mian_Fig_2aa


### perMAONVA
permanova_data1 = as.data.frame(permanova_mod$aov.tab)[1:10,]
permanova_data1$Label = c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                          "Year × Site × Origin", "Years × Species", "Site × Species")
permanova_data1$Label = factor(permanova_data1$Label, levels = rev(c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                                                     "Years × Species", "Site × Species", "Year × Site × Origin")))
permanova_data1$R2_per = sprintf("%.1f", permanova_data1$R2*100)
permanova_data1$R2 = round(permanova_data1$R2,2)
permanova_data1$`q-vaules`=p.adjust(permanova_data1$`Pr(>F)`, method = "BH")
permanova_data1$p_label = paste0("p = ", sprintf("%.3f", permanova_data1$`q-vaules`))
permanova_data1$p_shape = as.factor(ifelse(permanova_data1$`q-vaules` >= 0.05, 1, 0))


ggplot(permanova_data1, aes(x = Label, y = as.numeric(R2_per))) +  # 棒棒图的连线  
  geom_segment(aes(x = Label, xend = Label, y = 0, yend = as.numeric(R2_per)-0.015),                 
               linetype = "solid",size = 3,color = "#D3D3D3") + 
  geom_point(aes(shape = p_shape, fill = p_shape), color = "black", size = 10) + 
  geom_text(aes(label = R2_per),color = "black",size = 3) + 
  geom_text(aes(label = p_label),            
            hjust = ifelse(permanova_data1$R2_per >= 0, -0.4, 0.5),
            vjust = 0.2, angle = 0, color = "black",size = 3.5) + 
  labs(x = NULL, y = "Explained variance (%)", tag = "b") + 
  scale_shape_manual(values = c(22,22)) + 
  scale_fill_manual(values = c("#08ADB5","white")) + 
  theme_bw() + mytheme +
  theme(axis.text.x = element_text(size = 11, angle = 0),
        plot.margin = margin(l = 60, r = 20, t = 10, b = 10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0,30)) +
  coord_flip() -> mian_Fig_2b; mian_Fig_2b


############################## dispersion difference ###########################
rownames(as.matrix(Bray_dist_field_no)) %in% Field_group$Sample_ID
bdisper_field <- vegan::betadisper(Bray_dist_field_no, Field_group$Origin, type = "centroid", bias.adjust = FALSE)
## Test significance of dispersion difference
anova(bdisper_field)
TukeyHSD(bdisper_field)

Field_aov_result_all = NULL
Field_beita_data_all = NULL

Years = unique(Field_group$Years)
Site = unique(Field_group$Site)

## 分时间站点分别比较外来植物与本地植物种内组成相似性
for(i in Years){
  for(ii in Site) {
    ## Field survey
    select_group = subset(Field_group, Years == i & Site == ii)
    field_select_dist = as.matrix(Bray_dist_field_no)[select_group$Sample_ID, select_group$Sample_ID]
    select_group = select_group[rownames(field_select_dist),]
    bdisper_field <- betadisper(as.dist(field_select_dist), select_group$Origin, type = "centroid", bias.adjust = FALSE)
    ## Test significance of dispersion difference
    field_aov_mod = anova(bdisper_field)
    # TukeyHSD(bdisper_field)
    Field_aov_result = data.frame(Type = "Field", Years = i, Site = ii, F_value = field_aov_mod$`F value`[1], p_value = field_aov_mod$`Pr(>F)`[1])
    ##
    Field_beita_data = as.data.frame(bdisper_field$distances) 
    Field_beita_data$Sample_ID = rownames(Field_beita_data)
    colnames(Field_beita_data)[1] = "Field_beita"
    Field_beita_data = Field_beita_data %>% left_join(select_group[,c("Sample_ID", "Species")])
    Field_beita_data = Field_beita_data[,c("Sample_ID","Species","Field_beita")]
    ## 合并数据集
    Field_aov_result_all = rbind(Field_aov_result_all, Field_aov_result)
    Field_beita_data_all = rbind(Field_beita_data_all, Field_beita_data)
  }
}

###
Field_beita_data = Field_beita_data_all %>% left_join(Field_group[,c("Sample_ID","Species","Origin", "Years", "Site", "Latitude")])
head(Field_beita_data)
Beita_sum_data = Rmisc::summarySE(Field_beita_data, measurevar = "Field_beita", 
                                  groupvars = c("Years","Site","Origin"))

Beita_sum_data$Site = factor(Beita_sum_data$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Beita_sum_data$Origin = factor(Beita_sum_data$Origin, levels = c("Native","Exotic"))

# 创建灰色背景的数据框
background_data <- data.frame(
  Site = c("Guangzhou", "Changsha", "Zhengzhou"), 
  xmin = c(0.5, 2.5, 4.5),         
  xmax = c(1.5, 3.5, 5.5))

ggplot(Beita_sum_data, aes(x = Site, y = Field_beita, fill = Years, color = Years, shape = Origin)) + 
  geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), 
            fill = "gray95", alpha = 0.5, inherit.aes = FALSE) + 
  geom_point(aes(fill = Years),position = position_dodge(0.9), color = "black", size = 3) +
  scale_color_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  scale_fill_manual(values = c("#898EA1","#CF9742","#3A7C72")) +
  geom_errorbar(aes(ymin = Field_beita - se, ymax = Field_beita + se, color = Years), 
                width = 0, linewidth = 0.4, position = position_dodge(width = 0.9)) +
  geom_segment(aes(x = 1.9, xend = 2.1, y = 0.59, yend = 0.59), color = "black") + 
  annotate("text", x = 2, y = 0.60, label = expression(italic("p") ~ "= 0.035"), size = 3.5) +
  geom_segment(aes(x = 2.9, xend = 3.1, y = 0.62, yend = 0.62), color = "black") + 
  annotate("text", x = 3, y = 0.63, label = expression(italic("p") ~ "= 0.012"), size = 3.5) +
  geom_segment(aes(x = 3.2, xend = 3.4, y = 0.60, yend = 0.60), color = "black") + 
  annotate("text", x = 3.3, y = 0.61, label = expression(italic("p") ~ "= 0.026"), size = 3.5) +
  geom_segment(aes(x = 3.5, xend = 3.7, y = 0.62, yend = 0.62), color = "black") + 
  annotate("text", x = 3.6, y = 0.63, label = expression(italic("p") ~ "= 0.021"), size = 3.5) +
  geom_segment(aes(x = 6.2, xend = 6.4, y = 0.57, yend = 0.57), color = "black") + 
  annotate("text", x = 6.3, y = 0.58, label = expression(italic("p") ~ "= 0.001"), size = 3.5) +
  scale_shape_manual(values = c(21, 22)) +
  theme_bw() + mytheme + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) +   
  scale_y_continuous(limits =  c(0.34, 0.63)) + 
  labs(x = '', y = '') -> mian_Fig_2c; mian_Fig_2c


# 找到 p5 的 y 轴范围
y_limits <- c(min(Beita_sum_data$Field_beita - Beita_sum_data$se), 
              max(Beita_sum_data$Field_beita + Beita_sum_data$se))

library(Rmisc)  
data_all = Rmisc::summarySE(Beita_sum_data, groupvars = "Origin", measurevar = "Field_beita")
ggplot(data_all, aes(x = Origin, y = Field_beita, shape = Origin)) + 
  geom_point(fill = "black", color = "black", size = 3) +
  geom_errorbar(aes(ymin = Field_beita - se, ymax = Field_beita + se, group = Origin), 
                width = 0, linewidth = 0.4, color = "black") +
  geom_segment(aes(x = 1, xend = 2, y = 0.55, yend = 0.55)) + 
  annotate("text", x = 1.5, y = 0.56, label = expression(italic("p") ~ "= 0.062"), size = 3.5) +
  scale_shape_manual(values = c(21, 22)) +
  theme_bw() + mytheme + theme(legend.position = "none") +   
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        strip.placement = "outside",
        strip.text = element_text(size = 11)) + 
  scale_y_continuous(limits = y_limits) +  
  labs(x = '', y = 'Compositioanl dispersions ', tag = "c") -> mian_Fig_2cc; mian_Fig_2cc


(mian_Fig_2cc|mian_Fig_2c) + plot_layout(widths = c(0.2,0.8)) -> mian_Fig_2C; mian_Fig_2C

##################################
colnames(Field_beita_data_all)
beita_data_Centroid = merge(Field_beita_data_all[,c(1,3)], Field_group_scale, by = 'Sample_ID', all.x = TRUE)
colnames(beita_data_Centroid)
options(na.action ="na.fail")
colnames(Field_group)
pd_attributes_variable2 <- attributes(scale(beita_data_Centroid[c("Field_beita")]))


library(nlme)
multi_model = lmer(scale(Field_beita) ~ Site_pool + Soil_ph + Wcont + Soil_N + Tave + Prec + Funct_PC1 + Funct_PC2 + Phylo_PC1 + Phylo_PC2 +   
                     Funct_PC1:Soil_ph + Funct_PC1:Wcont + Funct_PC1:Soil_N + Funct_PC1:Tave + Funct_PC1:Prec + 
                     Funct_PC2:Soil_ph + Funct_PC2:Wcont + Funct_PC2:Soil_N + Funct_PC2:Tave + Funct_PC2:Prec + 
                     Phylo_PC1:Soil_ph + Phylo_PC1:Wcont + Phylo_PC1:Soil_N + Phylo_PC1:Tave + Phylo_PC1:Prec +
                     Phylo_PC2:Soil_ph + Phylo_PC2:Wcont + Phylo_PC2:Soil_N + Phylo_PC2:Tave + Phylo_PC2:Prec + 
                     (1|Years), REML=FALSE, data = beita_data_Centroid)

shapiro.test(residuals(multi_model))
car::Anova(multi_model)
summary(multi_model)
MuMIn::r.squaredGLMM(multi_model)

# Simplify the model - use Log-likelihood ration test backed up by AIC
drop1(multi_model, test = "Chi") # 
f.sbs1 <- update(multi_model, ~. -Soil_ph:Funct_PC1)
AIC(f.sbs1)
drop1(f.sbs1, test = "Chi")
f.sbs2 <- update(f.sbs1, ~. -Soil_ph:Phylo_PC1)
AIC(f.sbs2)
drop1(f.sbs2, test = "Chi")
f.sbs3 <- update(f.sbs2, ~. -Wcont:Funct_PC1)
AIC(f.sbs3)
drop1(f.sbs3, test = "Chi")
f.sbs4 <- update(f.sbs3, ~. -Tave:Phylo_PC2)
AIC(f.sbs4)
drop1(f.sbs4, test = "Chi")
f.sbs5 <- update(f.sbs4, ~. -Soil_ph:Phylo_PC2)
AIC(f.sbs5)
drop1(f.sbs5, test = "Chi")
f.sbs6 <- update(f.sbs5, ~. -Soil_N:Phylo_PC2)
AIC(f.sbs6)
drop1(f.sbs6, test = "Chi")
f.sbs7 <- update(f.sbs6, ~. -Prec:Funct_PC1)
AIC(f.sbs7)
drop1(f.sbs7, test = "Chi")
f.sbs8 <- update(f.sbs7, ~. -Wcont:Phylo_PC1)
drop1(f.sbs8, test = "Chi")
f.sbs9 <- update(f.sbs8, ~. -Prec:Funct_PC2)
drop1(f.sbs9, test = "Chi")
f.sbs10 <- update(f.sbs9, ~. -Tave:Funct_PC2)
drop1(f.sbs10, test = "Chi")
f.sbs11 <- update(f.sbs10, ~. -Wcont:Phylo_PC2)
drop1(f.sbs11, test = "Chi")
f.sbs12 <- update(f.sbs11, ~. -Soil_N:Phylo_PC1)
drop1(f.sbs12, test = "Chi")
f.sbs13 <- update(f.sbs12, ~. -Soil_N:Funct_PC2)
drop1(f.sbs13, test = "Chi")
f.sbs14 <- update(f.sbs13, ~. -Wcont:Funct_PC2)
drop1(f.sbs14, test = "Chi")
f.sbs15 <- update(f.sbs14, ~. -Wcont)
drop1(f.sbs15, test = "Chi")
f.sbs16 <- update(f.sbs15, ~. -Soil_ph:Funct_PC2)
drop1(f.sbs16, test = "Chi")
f.sbs17 <- update(f.sbs16, ~. -Funct_PC2)
drop1(f.sbs17, test = "Chi")
f.sbs18 <- update(f.sbs17, ~. -Prec:Phylo_PC2)
AIC(f.sbs18)
drop1(f.sbs18, test = "Chi")
f.sbs19 <- update(f.sbs18, ~. -Phylo_PC2)
AIC(f.sbs19)
drop1(f.sbs19, test = "Chi")
f.sbs20 <- update(f.sbs19, ~. -Soil_ph)
AIC(f.sbs20)
drop1(f.sbs20, test = "Chi")

library(flextable)
library(dplyr)
# needs the packages flextable and officer to be installed

my_theme_flextable=function (x) 
{
  if (!inherits(x, "flextable")) 
    stop("my_theme supports only flextable objects.")
  x <- border_remove(x)
  x <- hline_top(x, part = "header", border = officer::fp_border())
  x <- hline_bottom(x, part = "header", border = officer::fp_border())
  x <- hline_bottom(x, part = "body", border = officer::fp_border())
  x <- font(x,fontname = "Times")
  x <- font(x,part="header",fontname = "Times")
  x
}


library(MuMIn)
table_anova=data.frame(Step=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"),
                       Model=c("Full model","-Soil_ph:Funct_PC1","-Soil_ph:Phylo_PC1","-Wcont:Funct_PC1","-Tave:Phylo_PC2","-Soil_ph:Phylo_PC2","-Soil_N:Phylo_PC2","-Prec:Funct_PC1", "-Wcont:Phylo_PC1", "-Prec:Funct_PC2",
                               "-Tave:Funct_PC2", "-Wcont:Phylo_PC2", "-Soil_N:Phylo_PC1", "-Soil_N:Funct_PC2","-Wcont:Funct_PC2","-Wcont", "-Soil_ph:Funct_PC2","-Funct_PC2","-Prec:Phylo_PC2", "-Phylo_PC2","-Soil_ph"),
                       #n_par=c(anova(multi_model,f.sbs1)[1,1],anova(multi_model,f.sbs1)[2,1],anova(f.sbs1,f.sbs2)[2,1],anova(f.sbs2,f.sbs3)[2,1],anova(f.sbs3,f.sbs4)[2,1],anova(f.sbs4,f.sbs5)[2,1],anova(f.sbs5,f.sbs6)[2,1],anova(f.sbs6,f.sbs7)[2,1],anova(f.sbs7,f.sbs8)[2,1],anova(f.sbs8,f.sbs9)[2,1],anova(f.sbs9,f.sbs10)[2,1]),
                       Chisq=round(c(NA,anova(multi_model,f.sbs1)[2,6],anova(f.sbs1,f.sbs2)[2,6],anova(f.sbs2,f.sbs3)[2,6],anova(f.sbs3,f.sbs4)[2,6],anova(f.sbs4,f.sbs5)[2,6],anova(f.sbs5,f.sbs6)[2,6],anova(f.sbs6,f.sbs7)[2,6],anova(f.sbs7,f.sbs8)[2,6],anova(f.sbs8,f.sbs9)[2,6],anova(f.sbs9,f.sbs10)[2,6],
                                     anova(f.sbs10,f.sbs11)[2,6],anova(f.sbs11,f.sbs12)[2,6],anova(f.sbs12,f.sbs13)[2,6],anova(f.sbs13,f.sbs14)[2,6],anova(f.sbs14,f.sbs15)[2,6],anova(f.sbs15,f.sbs16)[2,6],anova(f.sbs16,f.sbs17)[2,6],anova(f.sbs17,f.sbs18)[2,6],anova(f.sbs18,f.sbs19)[2,6],anova(f.sbs19,f.sbs20)[2,6]),1),
                       df=round(c(NA,anova(multi_model,f.sbs1)[2,7],anova(f.sbs1,f.sbs2)[2,7],anova(f.sbs2,f.sbs3)[2,7],anova(f.sbs3,f.sbs4)[2,7],anova(f.sbs4,f.sbs5)[2,7],anova(f.sbs5,f.sbs6)[2,7],anova(f.sbs6,f.sbs7)[2,7],anova(f.sbs7,f.sbs8)[2,7],anova(f.sbs8,f.sbs9)[2,7],anova(f.sbs9,f.sbs10)[2,7],
                                  anova(f.sbs10,f.sbs11)[2,7],anova(f.sbs11,f.sbs12)[2,7],anova(f.sbs12,f.sbs13)[2,7],anova(f.sbs13,f.sbs14)[2,7],anova(f.sbs14,f.sbs15)[2,7],anova(f.sbs15,f.sbs16)[2,7],anova(f.sbs16,f.sbs17)[2,7],anova(f.sbs17,f.sbs18)[2,7],anova(f.sbs18,f.sbs19)[2,7],anova(f.sbs19,f.sbs20)[2,7]),1),
                       p_value=as.character(signif(c(NA,anova(multi_model,f.sbs1)[2,8],anova(f.sbs1,f.sbs2)[2,8],anova(f.sbs2,f.sbs3)[2,8],anova(f.sbs3,f.sbs4)[2,8],anova(f.sbs4,f.sbs5)[2,8],anova(f.sbs5,f.sbs6)[2,8],anova(f.sbs6,f.sbs7)[2,8],anova(f.sbs7,f.sbs8)[2,8],anova(f.sbs8,f.sbs9)[2,8],anova(f.sbs9,f.sbs10)[2,8],
                                                     anova(f.sbs10,f.sbs11)[2,8],anova(f.sbs11,f.sbs12)[2,8],anova(f.sbs12,f.sbs13)[2,8],anova(f.sbs13,f.sbs14)[2,8],anova(f.sbs14,f.sbs15)[2,8],anova(f.sbs15,f.sbs16)[2,8],anova(f.sbs16,f.sbs17)[2,8],anova(f.sbs17,f.sbs18)[2,8],anova(f.sbs18,f.sbs19)[2,8],anova(f.sbs19,f.sbs20)[2,8]),3)),
                       R2m=round(c(r.squaredGLMM(multi_model)[1],r.squaredGLMM(f.sbs1)[1],r.squaredGLMM(f.sbs2)[1],r.squaredGLMM(f.sbs3)[1],r.squaredGLMM(f.sbs4)[1],r.squaredGLMM(f.sbs5)[1],r.squaredGLMM(f.sbs6)[1],r.squaredGLMM(f.sbs7)[1],r.squaredGLMM(f.sbs8)[1],r.squaredGLMM(f.sbs9)[1],r.squaredGLMM(f.sbs10)[1],
                                   r.squaredGLMM(f.sbs11)[1],r.squaredGLMM(f.sbs12)[1],r.squaredGLMM(f.sbs13)[1],r.squaredGLMM(f.sbs14)[1],r.squaredGLMM(f.sbs15)[1],r.squaredGLMM(f.sbs16)[1],r.squaredGLMM(f.sbs17)[1],r.squaredGLMM(f.sbs18)[1],r.squaredGLMM(f.sbs19)[1],r.squaredGLMM(f.sbs20)[1]),2),
                       R2c=round(c(r.squaredGLMM(multi_model)[2],r.squaredGLMM(f.sbs1)[2],r.squaredGLMM(f.sbs2)[2],r.squaredGLMM(f.sbs3)[2],r.squaredGLMM(f.sbs4)[2],r.squaredGLMM(f.sbs5)[2],r.squaredGLMM(f.sbs6)[2],r.squaredGLMM(f.sbs7)[2],r.squaredGLMM(f.sbs8)[2],r.squaredGLMM(f.sbs9)[2],r.squaredGLMM(f.sbs10)[2],
                                   r.squaredGLMM(f.sbs11)[2],r.squaredGLMM(f.sbs12)[2],r.squaredGLMM(f.sbs13)[2],r.squaredGLMM(f.sbs14)[2],r.squaredGLMM(f.sbs15)[2],r.squaredGLMM(f.sbs16)[2],r.squaredGLMM(f.sbs17)[2],r.squaredGLMM(f.sbs18)[2],r.squaredGLMM(f.sbs19)[2],r.squaredGLMM(f.sbs20)[2]),2),
                       AIC=round(c(AIC(multi_model),AIC(f.sbs1),AIC(f.sbs2),AIC(f.sbs3),AIC(f.sbs4),AIC(f.sbs5),AIC(f.sbs6),AIC(f.sbs7),AIC(f.sbs8),AIC(f.sbs9),AIC(f.sbs10),
                                   AIC(f.sbs11),AIC(f.sbs12),AIC(f.sbs13),AIC(f.sbs14),AIC(f.sbs15),AIC(f.sbs16),AIC(f.sbs17),AIC(f.sbs18),AIC(f.sbs19),AIC(f.sbs20)),1),
                       AIC_weight=round(MuMIn::Weights(AIC(multi_model,f.sbs1,f.sbs2,f.sbs3,f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10,
                                                           f.sbs11,f.sbs12,f.sbs13,f.sbs14,f.sbs15,f.sbs16,f.sbs17,f.sbs18,f.sbs19,f.sbs20)), 3)) %>% as.data.frame() %>%
  mutate(delta_AIC=round(AIC-min(AIC),1)) %>%
  dplyr::select(Step:AIC,delta_AIC,AIC_weight) %>%
  flextable() %>%
  bold(i=~delta_AIC==0) %>%
  set_header_labels(Chisq="χ²",p_value="p-value",R2c="R²c",R2m="R²m",delta_AIC="∆AIC",AIC_weight="AIC weight" ) %>%
  my_theme_flextable() %>%
  autofit()
table_anova


library(lmerTest)
f.sbs.fin <- update(f.sbs19, ~., REML= TRUE) # Update to REML to extract estimates
summary(f.sbs.fin)
anova(f.sbs.fin)
ranova(f.sbs.fin)
MuMIn::r.squaredGLMM(f.sbs.fin)
effectsize::effectsize(f.sbs.fin)


## Table (Fixed trem)
Table_mod = as.data.frame(anova(f.sbs.fin))
Table_mod$p = round(Table_mod$`Pr(>F)`, 3)
Table_mod$p.adj=round(p.adjust(Table_mod$`Pr(>F)`, method = "BH"), 3)
#Table_mod$p.adj=ifelse(Table_mod$p.adj < 0.001, "< 0.001", Table_mod$p.adj)
Table_mod$`Fixed factors` = c("Site pool",
                              "Soil pH", "Soil N", 
                              "MAT", "MAP", 
                              "FunctPC1", "PhyloPC1",
                              "Soil N × FunctPC1", 
                              "MAT × FunctPC1", 
                              "MAT × Phylo_PC1", 
                              "MAP × Phylo_PC1")

Table_mod$df = paste0(Table_mod$NumDF, ", ", round(Table_mod$DenDF,0))
Table_mod$`F` = round(Table_mod$`F value`, 2)
Table_mod$VIF = round(vif(f.sbs.fin), 2)
rownames(Table_mod) = NULL
print(Table_mod[,c("Fixed factors","df","F","p","p.adj","VIF")])
Table_mod_plot = Table_mod[,c("Fixed factors","df","F","p","p.adj","VIF")] %>%
  flextable() %>% 
  bold(j = "p.adj",  # 指定要加粗的列
       i = ~ p.adj <= 0.05,  # 加粗条件
       part = "body") %>%  # 应用于表格主体部分
  my_theme_flextable() %>% autofit()
print(Table_mod_plot)


### p.adjust
tables2<- as.data.frame(coef(summary(f.sbs.fin)))
tables2=tables2[-which(rownames(tables2) == "(Intercept)"),]
tables2$`p-value` = round(tables2$`Pr(>|t|)`, 3)
tables2$`q-vaules`=p.adjust(tables2$`Pr(>|t|)`, method = "BH")
tables2$`q-vaules` = round(tables2$`q-vaules`, 3)
tables2$Parameter = rownames(tables2)

################################################################################
library(dplyr)
library(ggplot2)
MegaModelSummary  = as.data.frame(coef(summary(f.sbs.fin)))[-1,]
MegaModelSummary$Parameter = rownames(MegaModelSummary)
MegaModelSummary$Parameter2 = c("Site pool","Soil pH", "Soil N", "MAT", "MAP", "FunctPC1", "PhyloPC1",
                                "Soil N × FunctPC1", "MAT × FunctPC1", "MAT × Phylo_PC1", "MAP × Phylo_PC1")
MegaModelSummary$Parameter2 = factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Soil pH", "Soil N", "MAT", "MAP", "FunctPC1", "PhyloPC1",
                                                                                 "Soil N × FunctPC1", "MAT × FunctPC1", "MAT × Phylo_PC1", "MAP × Phylo_PC1")))
MegaModelSummary = MegaModelSummary %>% left_join(tables2[,c("Parameter","q-vaules")], by = "Parameter")
MegaModelSummary$p_label = paste0("p = ", sprintf("%.3f", MegaModelSummary$`q-vaules`))
MegaModelSummary$sig2 <- as.factor(ifelse(MegaModelSummary$`q-vaules` > 0.05, 2, 1))

MegaModelSummary$Group = c("Site pool", "Soil properties", "Soil properties", "Climate", "Climate",
                           "Plant attributes", "Plant attributes","Interaction", "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group = factor(MegaModelSummary$Group, levels = c("Site pool", "Soil properties", "Climate", "Plant attributes", "Interaction"))


ggplot(MegaModelSummary, aes(x = Parameter2, y = Estimate, fill = Group))+
  geom_errorbar(aes(ymin=Estimate-1.96*`Std. Error`, ymax=Estimate+1.96*`Std. Error`), width=0, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 11.2), color = "black", linetype = "dashed") + 
  geom_text(aes(y = Estimate+1.96*`Std. Error`, 
                label = ifelse(is.na(p_label), "", paste("italic(p)==", `q-vaules`))),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.5) + 
  annotate("text", x = 12, y = 0.01,
           label = expression("Best model: " * italic(R)^2 * "m = 0.15, " * italic(R)^2 * "c = 0.16"),
           colour = "black", size = 4) + 
  labs(x = '', y = 'Standard regression coefficients', color = '', tag = "d") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#45A0A3","#88D8BF", "#F3DBC1","#9796C8","#473C8B")) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x =  element_text(color = "black", size = 13),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_shape_manual(values = c(16,21)) -> mian_Fig_2d; mian_Fig_2d

library(forcats)
MegaModelSummary_var = MegaModelSummary %>%
  mutate(Group=fct_relevel(Group,
                           (c("Site pool","Soil properties",
                              "Climate","Plant attributes",
                              "Interaction")))) %>% 
  group_by(Group) %>%
  summarise(sum_value=sum(abs(Estimate))) %>% 
  mutate(new_col=sum_value/sum(sum_value)) 

MegaModelSummary_var$Group = factor(MegaModelSummary_var$Group, levels = c("Site pool", "Soil properties", "Climate", "Plant attributes", "Interaction"))

ggplot()+
  geom_bar(data = MegaModelSummary_var, aes(x = "", y = new_col, fill = Group), 
           stat = "identity", width = 0.5, color = "black")+
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("#45A0A3","#88D8BF", "#F3DBC1","#9796C8","#473C8B")) +
  #scale_fill_manual(values = c("#45A0A3","#88D8BF", "#F3DBC1","#E68C51","#D75C4D")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA))+
  guides(fill = guide_legend(title = "Province",ncol = 3, nrow = 2, override.aes = list(shape = 22, size=0.2))) +
  labs(y = "Relative effect of estimates (%)") -> mian_Fig_2dd; mian_Fig_2dd

library(patchwork)
(mian_Fig_2d+mian_Fig_2dd) + plot_layout(widths = c(0.7,0.3)) -> mian_Fig_2D; mian_Fig_2D


################################################################################
library(ggeffects)
predmw.m1.rsc2=ggeffect(f.sbs.fin, terms = c("Tave","Funct_PC1"))
plot(predmw.m1.rsc2)

eff_mod_data <- data.frame(predmw.m1.rsc2)
colnames(eff_mod_data)[1] = "Tave"
colnames(eff_mod_data)[2] = "Field_beita"
colnames(eff_mod_data)[6] = "Funct_PC1"
# backtransform pd_attributes_variable2
eff_mod_data["Tave"] <- pd_attributes_variable$`scaled:center`["Tave"] + pd_attributes_variable$`scaled:scale`["Tave"]*eff_mod_data["Tave"]
eff_mod_data["Field_beita"] <- pd_attributes_variable2$`scaled:center`["Field_beita"] + pd_attributes_variable2$`scaled:scale`["Field_beita"]*eff_mod_data["Field_beita"]

################################################################################
library(ggeffects)
predmw.m1.rsc2=ggeffect(f.sbs.fin, terms = c("Funct_PC1","Tave"))
plot(predmw.m1.rsc2)

eff_mod_data <- data.frame(predmw.m1.rsc2)
colnames(eff_mod_data)[1] = "Funct_PC1"
colnames(eff_mod_data)[2] = "Field_beita"
colnames(eff_mod_data)[6] = "Tave"
# backtransform pd_attributes_variable2
eff_mod_data["Funct_PC1"] <- pd_attributes_variable$`scaled:center`["Funct_PC1"] + pd_attributes_variable$`scaled:scale`["Funct_PC1"]*eff_mod_data["Funct_PC1"]
eff_mod_data["Field_beita"] <- pd_attributes_variable2$`scaled:center`["Field_beita"] + pd_attributes_variable2$`scaled:scale`["Field_beita"]*eff_mod_data["Field_beita"]

eff_mod_data$Tave <- ifelse(eff_mod_data$Tave == "-1", "- 1 SD", 
                            ifelse(eff_mod_data$Tave == "0", "Mean", "+ 1 SD"))
eff_mod_data$Tave <- factor(eff_mod_data$Tave, levels = c("- 1 SD", "Mean", "+ 1 SD"))

ggplot()+
  #geom_point(data = lme_mod_data2, mapping = aes(x = Tave_row, y = Effect_size), pch = 21) + 
  geom_line(data = eff_mod_data, mapping = aes(Funct_PC1,Field_beita,color=factor(Tave)), size=1.5)+
  #geom_ribbon(data = eff_mod_data, mapping = aes(x = Funct_PC1, ymin=conf.low,ymax=conf.high,fill=factor(Funct_PC1)),alpha=0.3, colour= NA)+
  #scale_fill_manual(values=c("#FEEC00","#C0B667", "#737891", "#425C8E", "#001B56"),guide=F)+
  #scale_color_manual(values=c("#FEEC00","#C0B667", "#737891", "#425C8E", "#001B56"),breaks=waiver(), name="Phylo Di (log10)")+
  theme_bw() + mytheme + 
  guides(color = guide_legend(override.aes = list(fill = NA)))+
  scale_fill_manual(values=rev(c("#324554", "#B89088", "#E9D1CD")),guide=F)+
  scale_color_manual(values=rev(c("#324554", "#B89088", "#E9D1CD")),breaks=waiver(), name="Mean annual temperature")+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  #annotate("text", label = expression(italic(p) == 0.024), x = 15, y = 0.25, size = 4) + 
  #scale_y_continuous(limits = c(0.47, 0.59),breaks = seq(0.47, 0.59, by = 0.04)) + 
  labs(x = "The PC1 of plant traits", 
       y = "Compositional dispersions\n based on Bray-Curtis distance", tag = "e", color = "Mean annual temperature") +
  annotate("segment", x = -2, xend = -4, y = 0.52, 
           yend = 0.52, colour = "black", size = 0.5, arrow = arrow(angle = 15, 
                                                                    length = unit(0.5,  "cm"))) + 
  annotate("segment", x = 0, xend = 2, y = 0.52, 
           yend = 0.52, colour = "black", size = 0.5, arrow = arrow(angle = 15, 
                                                                    length = unit(0.5,  "cm"))) +
  theme(legend.position = c(0.65,0.18),
        legend.key = element_blank(),
        legend.title = element_text(size=11),
        legend.text= element_text(size=10),
        legend.background = element_blank())-> mian_Fig_2E; mian_Fig_2E


### 13.20 x 5.39
(mian_Fig_2D|mian_Fig_2E) + plot_layout(widths = c(0.62,0.38))






