library(phytools)
library(ggtree)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(sf)
library(ggplot2)
library(ggspatial)

###
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="NA"),
                   axis.line.y=element_line(size=0.5, colour="NA"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=11),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=9),
                   legend.title = element_text(size=10),
                   plot.tag = element_text(size = 14, face = "bold"))

## 采样点绘制
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") 
class(china_map)

### sample site information
datasel = data.frame(Site = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"),
                     Latitude = c(23.1, 25.2, 27.9, 30.5, 34.6, 36.2),
                     Longitude = c(113.2, 110.2, 112.9, 114.3, 113.6, 117.1))


ggplot(china_map)+
  geom_sf(data = china_map,fill="grey95",size=1) + 
  xlim(105,122)+ ylim(18,42)+ 
  ggnewscale::new_scale_fill() + 
  #geom_text(data=datasel,aes(x=lon, y=lat-0.6 ,label=province),size=3,colour="black")+
  geom_point(data = datasel, aes(x=Longitude ,y=Latitude),
             size=4,alpha=1, shape = 21, color = "black", fill = "grey50") + 
  main_theme +
  annotation_scale(location = "br", style = "ticks",line_width = 0.1,pad_y = unit(0.5, "cm"), text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"),width = unit(1, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),style = north_arrow_fancy_orienteering) +
  #coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=104")+
  guides(fill = guide_legend(title = "Site",ncol = 1, nrow = 14, override.aes = list(shape = 22, size=5))) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill='#FFFFFF', colour='black'),
        axis.line.x = element_line(size=0.5, colour="black"),
        axis.line.x.top =element_line(size=0.5, colour="black"),
        axis.line.y.right = element_line(size=0.5, colour="black"),
        axis.line.y=element_line(size=0.5, colour="black"),
        legend.position = "right",
        panel.grid.major = element_line(color = "white",size=0.2))+
  labs(x='', y='') +
  ggrepel::geom_text_repel(mapping = aes(x=Longitude,y=Latitude,label=Site), data = datasel,size = 3.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25)) -> p1; p1

## 采样样本数统计
# Soil sample grouping information
Field_group = read.xlsx("field_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

sample_num_sum = Field_group %>% group_by(Years, Site, Origin) %>%
  summarise(sample_num = n())

df_wide <- sample_num_sum %>%
  tidyr::pivot_wider(names_from = Origin, values_from = sample_num)
df_wide$Total = df_wide$Exotic + df_wide$Native
df_wide$Site = factor(df_wide$Site, levels = rev(c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an")))


df_long <- df_wide %>%
  tidyr::pivot_longer(cols = c(Exotic, Native, Total),  # 选择需要转换的列
                      names_to = "Group",       # 新列的名称
                      values_to = "sample_num")  # 新值的列名

##
df_long$Years = as.factor(df_long$Years)
df_long$Group = factor(df_long$Group, levels = c("Native","Exotic","Total"))

ggplot(df_long, aes(Group, Years)) +
  geom_tile(color = "black", fill = NA) +
  #scale_fill_gradient2(low = ("#A67C2A"),mid = "white",high = muted("#001260"), midpoint = 0, limits=c(-0.68,0.85)) + # muted
  geom_text(aes(label = sample_num), size = 2.8, color = 'black') + 
  ggforce::facet_col(Site~., space = "free") +
  labs(x = NULL, y = NULL) + 
  theme(strip.text.x = element_blank(), 
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, angle = 35, vjust = 1, hjust = 1, color = "black"),
        panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = 'transparent')) -> p2


library(patchwork)
(p1|p2) + plot_layout(widths = c(0.9,0.18))


################################################################################
plant_tree = read.tree("IQ_tree_plant_2025.NEWICK")
plot(plant_tree)
## 分组信息
traits_mean = read.xlsx("Greenhouse_data_group.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
colnames(traits_mean)

tree_df <- tibble::as_tibble(plant_tree) 
tree_df$Species = tree_df$label
tree_df = tree_df %>% left_join(traits_mean[,c( "Species","Origin","Family","Genus")])

Family_color = c("Acanthaceae"="#332500","Amaranthaceae"="#542D20","Asteraceae" ="#994240","Caryophyllaceae" = "#D75C4D", "Cyperaceae" = "#E68C51",
                 "Euphorbiaceae" ="#F59D52", "Fabaceae" = "#EFBA55","Lamiaceae" ="#FCD170","Malvaceae"="#FEE1B0", "Onagraceae" = "#C5E0D0",
                 "Phytolaccaceae" = "#ABDCE0","Poaceae" ="#7FC5D8","Polygonaceae"="#73BCD5", "Solanaceae" = "#528FAC", "Urticaceae" = "#376694", "Verbenaceae" = "#1F466F", "NA" = "grey40")

tree_df$Origin = factor(tree_df$Origin, levels = c("Native", "Exotic"))
#to_drop<-c("Amborella_trichopoda","")
#plant_tree <- drop.tip(as.phylo(plant_tree), to_drop) 
ggtree(plant_tree,layout = "fan", branch.length="none",ladderize = F,size = 0.8,             
       open.angle =180) %<+% tree_df +
  geom_tiplab(aes(label = sub("_", " ", label)), size = 3, offset=0.3, fontface = "italic") +
  geom_tree(aes(color = Family), size = 0.5) + 
  scale_color_manual(values = Family_color) +
  ggnewscale::new_scale_color() + 
  geom_tippoint(aes(color = Origin, shape = Origin), size = 1.5) +
  #scale_color_manual(values = c("#2F9F84","#E8B66C")) + 
  scale_color_manual(values = c("#58C5BF","#FF8942")) +
  scale_shape_manual(values = c(16,15)) +
  theme(legend.position = "none") 


