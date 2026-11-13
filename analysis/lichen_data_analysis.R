                                     ### **************** ###
                                ### PLOTS INCLUDED  MANUSCRIPT ###
                                    ### **************** ###

############################################################################################
#####      1.     CLUSTERING sampling locations based on 6 bioclimatic variables       #####
############################################################################################

############################################################################################
#####      2.                    Plotting Taxonomic Composition                        #####
############################################################################################

############################################################################################
#####      3.    WEATHER and CLIMATE relationship with proteomic composition           #####
############################################################################################

############################################################################################
#####      4.       Plotting PERMANOVA results OF WEATHER and CLIMATE                  #####
############################################################################################

############################################################################################
#####      5.     VOLCANO PLot of DEPs and other NS proteins                           #####
############################################################################################



                                ### **************** ###
                             ###   MAY BE SUPPLEMENTARY  ###
                               ### **************** ###


############################################################################################
#####      6.          Differentially Expressed Proteins Only -  bar plot              #####
############################################################################################

############################################################################################
#####      7.        oxido-reductase activity proteins - GO Plot                       #####
############################################################################################

############################################################################################
#####      8.       Plotting level-1 molecular function GO terms of DEPs               #####  
#####                              DEPs only                                           #####
############################################################################################

############################################################################################
#####      9.                    ALGAL AND FUNGAL genotypes                            #####
############################################################################################


###### important libraries #####
library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
# Unused dependency omitted: reshape2
library(ggrepel)
library(ggpp)    
library(phyloseq)
library(fantaxtic)
library(factoextra)
# Unused dependency omitted: ggnested
library(knitr)
# Unused dependency omitted: ggforce
library(GOplot)
library(ape)
library(stringr)
                                
                                      ### **************** ###
                                ### PLOTS INCLUDED  MANUSCRIPT ###
                                     ### **************** ###

############################################################################################
#####   1.       CLUSTERING sampling locations based on 6 bioclimatic variables       #####
############################################################################################

# import data set and subset
data_env_prot <- read.csv("data_env_prot.csv", fileEncoding = "latin1")
env <- data_env_prot[,15:20] 
env <- as.data.frame(scale(env))
env=env/50                # dividing all the values of the cells by a specified number, e.g., by 50 in this case, 
                          # helps in shortening the dendrogram height
rownames(env) <- data_env_prot$sites_names

# Compute hierarchical k-means clustering
res.hk <- hkmeans(env, 4, hc.metric = "euclidean")  # define number of clusters

# visualize "rectangular", "dendrograms 
treePlot <- fviz_dend(res.hk,horiz = T,cex = 0.5, lwd = 0.5,k = 4,
                      k_colors = c("darkblue" ,"darkred", "black","darkgreen"),
                      color_labels_by_k = T,show_labels = T,repel = T,rect = F,
                      rect_border = "black",rect_lty = 1,rect_fill = T,
                      lower_rect = F,type = "rectangle",pylo_layout = "layout_as_tree",
                      labels_track_height = T,ggtheme = theme_void())


options(repr.plot.width=8, repr.plot.height=5)
plot(treePlot)        # visualize plot


############################################################################################
#####      2.                    Plotting Taxonomic Composition                       #####
############################################################################################

# https://github.com/gmteunisse/fantaxtic

## The package, fantaxtic,  was originally developed for a phyloseq object, however, it could also be used for 
## any tabular data set if a phyloseq object would be easily created using the data set. To create a phyloseq object
## we need the taxonomy information (TAXONOMY_data), the abundance information (ASV_data) and climatic region data (META_data)


# Creating a phyloseq object
# import all data sets
tax_abundance_data  <- read.csv("tax_abundance_data.csv", header=T, row.names=1, fileEncoding="latin1", check.names = FALSE)  ### fake ASVs 
tax_class_data  <- read.csv("tax_class_data.csv", header=T, row.names=1, fileEncoding="latin1", check.names = FALSE)  ### taxonomic abundance
meta_data <- read.csv("meta_data.csv", header=T, fileEncoding="latin1", check.names = FALSE) ### fake META 

## With phyloseq format
tax_abundance_data =   otu_table(data.frame(tax_abundance_data), taxa_are_rows = TRUE)
tax_class_data =   tax_table(as.matrix(tax_class_data))
meta_data    =   sample_data(data.frame(meta_data, row.names = meta_data$Sites))

## Make one phyloseq object, which contains ASVs and TAX objects
phyloseq_object <- phyloseq(tax_abundance_data, tax_class_data)

## Now create a random phylogenetic tree with the ape package, and add it to your phyloseq object created 
random_tree = rtree(ntaxa(phyloseq_object), rooted=TRUE, tip.label=taxa_names(phyloseq_object))


## Create the final phyloseq object which contains the meta data and the random tree above
phyloseq_object = merge_phyloseq(phyloseq_object, meta_data, random_tree)


# Clean data set (phyloseq)
phyloseq_Cleaned = subset_taxa(phyloseq_object, !is.na(Phylum) & !is.na(Family)) 

# Melt for ordering purposes
phyloseq_CleanedM=psmelt(phyloseq_Cleaned)

# Order taxa
phyloseq_CleanedM$Phylum=factor(phyloseq_CleanedM$Phylum, 
                         levels = c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                    "Chloroflexi","Firmicutes","Planctomycetes",
                                    "Proteobacteria","Tenericutes","Verrucomicrobia",
                                    "Other"))

# Order climatic regions
sample_data(phyloseq_Cleaned)$Climatic_Region<- factor(sample_data(phyloseq_Cleaned)$Climatic_Region, 
                                           levels = c("SAL","ALP","HM","SCL"))


## arrange which taxa to display (phylum and family in this case) 
top_nested <- nested_top_taxa(phyloseq_Cleaned,top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 9, n_nested_taxa = 5, 
                              nested_merged_label = "Other <tax>")

top_nested$top_taxa %>%
  mutate(top_abundance = round(top_abundance, 3),
         nested_abundance = round(nested_abundance, 3)) %>%
  kable(format = "markdown")


## visualize nested barplot with ggplot2 default colors of the package
plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Family",
                nested_merged_label = "Other <tax>",
                legend_title = "Phylum & Top 5 Family") +
  facet_grid(~Climatic_Region, scales = "free", space='free') +
  scale_x_discrete(expand = c(0, 0.5))+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(), 
        panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1, face="bold"),
        axis.text.y = element_text(angle=0, hjust=0.5, face="bold"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.2, "cm"),
        legend.title.align = 0)+
  guides(fill=guide_legend(ncol=1))


## visualize the plot with custom colors for phylums 
plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Family",
                nested_merged_label = "Other",
                palette = c(Acidobacteria = "lightblue", 
                            Actinobacteria = "red",
                            Firmicutes = "blue", 
                            Planctomycetes = "yellow",
                            Proteobacteria = "green"),
                merged_clr = "gray",
                legend_title = "Top 5 Phylum & Family")+
  facet_grid(~Climatic_Region, scales = "free", space='free') +
  scale_x_discrete(expand = c(0, 0.5))+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1, face="bold"),
        axis.text.y = element_text(angle=0, hjust=0.5, face="bold"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.title.align = 0)+
  guides(fill=guide_legend(ncol=1))



# Order samples by the total abundance of Proteobacteria
sample_order <- psmelt(top_nested$ps_obj) %>%
  data.frame() %>%
  
  # Calculate relative abundances
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  
  # Sort by taxon of interest
  filter(Phylum == "Proteobacteria") %>%
  group_by(Sample) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(Abundance) %>% 
  
  # Extract the sample order
  pull(Sample) %>%
  as.character()

# visualize plot which is ordered according to abundance of Proteobacteria
plot_nested_bar(top_nested$ps_obj,
                "Phylum",
                "Family",
                sample_order = sample_order,
                nested_merged_label = "Other",
                legend_title = "Phylum & Top 5 Family")+
  facet_grid(~Climatic_Region, scales = "free", space='free') +
  scale_x_discrete(expand = c(0, 0.5))+
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1, face="bold"),
        axis.text.y = element_text(angle=0, hjust=0.5, face="bold"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.title.align = 0)+
  guides(fill=guide_legend(ncol=1))

############################################################################################
#####      3.    WEATHER and CLIMATE relationship with proteomic composition           #####
############################################################################################

data_env_prot<-read.csv("data_env_prot.csv", fileEncoding = "latin1")
prot <- data_env_prot[,c(1,21:237)]
rownames(prot)<-prot$sites_names
prot=prot[,-1]

env_climate <- data_env_prot[,c(1, 15:20)]
rownames(env_climate)<-env_climate$sites_names
env_climate<-env_climate[,-1]
env_climate<-scale(env_climate)

env_weather <- data_env_prot[,c(1, 3:6)]
rownames(env_weather)<-env_weather$sites_names
env_weather<-env_weather[,-1]
env_weather<-scale(env_weather)

###  convert prot to a matrix ###
m_prot <- as.matrix(prot)

# nmds code
set.seed(123)
nmds <- metaMDS(m_prot, distance = "bray")

# run the envfit function with the climate and the local weather data frame, env_climate and env_local
en_climate <- envfit(nmds, env_climate, permutations = 999, na.rm = TRUE)
en_weather <-   envfit(nmds, env_weather,   permutations = 999, na.rm = TRUE)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores <- as.data.frame(scores(nmds)$sites)

# add 'season' column as before
data.scores$climatic_groups <- data_env_prot$climatic_groups

### extract information - new addition remove multiplication
en_climate_coord_cont <- as.data.frame(scores(en_climate, "vectors")) 
en_weather_coord_cont   <- as.data.frame(scores(en_weather,   "vectors")) 

data.scores$sites <- rep(c("Rhn", "SRa", "Drß", "Jhs",
                           "Kln", "Lme", "Rbd", "RSv",
                           "Sån", "Sns", "Vga", "Vbg", 
                           "Klp", "Tbg"), each=3)

### NMDS plot - climate and weather combined
climate_and_weather_NMDS <- ggplot() + geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2,
                                                                     color = data_env_prot$climatic_groups),size =5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color="gray40") + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color="gray40")+
  geom_segment(data = en_climate_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(angle = 22.5,length = unit(0.2,"cm"),type = "closed"),linetype = 1, linewidth = 0.6,colour = "gray40") +
  geom_segment(data = en_weather_coord_cont,aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(angle=22.5,length = unit(0.2,"cm"),type = "closed"),linetype=1, linewidth=0.6,colour = "skyblue3")+
  guides(color = guide_legend(title = "Climatic Regions", color = "black", hjust = 1)) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = sites, color = data_env_prot$climatic_groups),data = data.scores, size = 2, show.legend=FALSE) +
  geom_text_repel(data = en_climate_coord_cont, aes(x = NMDS1, y = NMDS2,
                                                    label = row.names(en_climate_coord_cont)), 
                  position = position_nudge_center(x = 0.001, center_x = 0,
                                                   y = 0.001, center_y = 0),
                  hjust = 0.5,size = 4, colour = "gray40",fontface = "bold") +
  geom_text_repel(data = en_weather_coord_cont,
                  aes(x = NMDS1, y = NMDS2,
                      label = row.names(en_weather_coord_cont)), 
                  position = position_nudge_center(x = 0.001, center_x = 0,
                                                   y = 0.001, center_y = 0),
                  hjust = 0.5,size = 4, colour = "skyblue",
                  fontface = "bold") +
  theme_minimal() +
  theme (axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.title.x = element_text(size = 12, hjust = 0.5),
         axis.title.y = element_text(size = 12, hjust = 0.5),
         plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
         panel.background = element_rect(fill = "white"),
         legend.text = element_text(colour = "black", size = 10, hjust=0),
         legend.title = element_text(colour = "Black", size = 11,face = "bold", hjust = 0),
         legend.position = c(0.90, 0.23),
         panel.border = element_rect(colour = "gray",fill = NA, linewidth = 0.5),
         strip.text = element_text(size = 18, face = "italic", hjust = 0.5),
         axis.text = element_text(size = 12),
         legend.key.size = unit(0.5, "cm"),
         legend.spacing.y = unit(0.5, "cm"),
         legend.background = element_rect(linetype = 1, linewidth = 0.5, colour = "lightgray"))+
  scale_colour_manual(values = c("darkblue" ,"black", "darkred", "darkgreen")) +
  scale_x_continuous(limits = c(-0.8, 0.8), expand = c(0.0, 0.0)) +
  scale_y_continuous(limits = c(-0.4, 0.4), expand = c(0.0, 0.0))


options(repr.plot.width=8, repr.plot.height=6)

plot(climate_and_weather_NMDS)


############################################################################################
#####       4.       Plotting PERMANOVA results OF WEATHER and CLIMATE                #####
############################################################################################

##### import data set #####
data_permanova_output<-read.csv("data_permanova_output.csv")  # data set with r squared values extracted from permanova analysis

# custom reorder the environmental variables
data_permanova_output$variable <- factor(data_permanova_output$variable, 
                      levels = c("Residual","Temperature Seasonality","Annual Mean Temperature",
                                "Mean Diurnal Range (Temperature)","Precipitation of Driest Quarter" ,
                                "Annual Precipitation","Precipitation Seasonality",
                                "Humidity_2W","Temp_2W","SampleDate_Julian","Precipcover_2W"))

# Plot the results - relative frequency bar chart
permPlot=ggplot(data_permanova_output, aes(x = catagory, y = rSquare, fill = variable)) +
  scale_fill_manual(values=c("#b3b3b3", "#dc1c13","#ea4c46", "#f07470", 
                             "#0000ff", "#4949ff", "#7879ff",
                             "#388688", "#7ebe99", "#a2d898","#cff39c"))+
  geom_bar(stat = 'identity') +
  labs(title = '', x = '', y = '', fill = 'Variables')+
  theme_bw()

plot(permPlot)

############################################################################################
#####         5.     VOLCANO PLot of DEPs and other NS proteins                        #####
############################################################################################

# DEP: Differentially expressed proteins
# NS: Not Significant

# import data set of the limma output
all_protein_limma_output=read.csv("all_protein_limma_output.csv", fileEncoding = "latin1")

# subset import proteins for labeling
data_label=  subset(all_protein_limma_output, label=="yes")  # for labeling 


# visualize plot with text labels in transparent box
volcanoPlot_boxLabel=ggplot(data=all_protein_limma_output, aes(x=logFC, y=B_statistics, size=neg_log10_p_val,col=Activity)) +
  geom_point(alpha=0.4) + theme_minimal() +
  geom_label_repel(aes(x=logFC, y=B_statistics, label=protein), 
                   size=2, min.segment.length = 0, seed = 42, alpha=0.5,box.padding = 0.5,
                   max.overlaps = Inf,arrow = arrow(length = unit(0.010, "npc")),
                   nudge_x = .15,nudge_y = .5,data = data_label, show.legend=F) +
  scale_color_manual(values=c("blue",  "red", "gray40")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5, color="gray40") + 
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.5, color="gray40")+
  theme (axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.title.x = element_text(size = 12, hjust = 0.5),
         axis.title.y = element_text(size = 12, hjust = 0.5),
         plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
         panel.background = element_rect(fill = "white"),
         legend.text = element_text(colour = "black", size = 10, hjust=0),
         legend.title = element_text(colour = "Black", size = 11,face = "bold", hjust = 0),
         legend.position = c(0.93, 0.19),
         panel.border = element_rect(colour = "gray",fill = NA, linewidth = 0.5),
         strip.text = element_text(size = 18, face = "italic", hjust = 0.5),
         axis.text = element_text(size = 12),
         legend.key.size = unit(0.5, "cm"),
         legend.spacing.y = unit(0.5, "cm"),
         legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgray"))

plot(volcanoPlot_boxLabel)

# visualize plot with text labels with out box
volcanoPlot_NO_boxLabel=ggplot(data=all_protein_limma_output, aes(x=logFC, y=B_statistics, size=neg_log10_p_val,col=Activity)) +
  geom_point(alpha=0.4) +theme_minimal() +
  geom_text_repel(aes(x=logFC, y=B_statistics, label=protein),
                  size = 2,min.segment.length = 0,seed = 42,box.padding = 0.5,
                  max.overlaps = Inf,arrow = arrow(length = unit(0.010, "npc")),
                  nudge_x = .15,nudge_y = .5,data = data_label, show.legend=F) +
  scale_color_manual(values=c("blue",  "red", "gray40")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5, color="gray40") + 
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.5, color="gray40")+
  theme (axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.title.x = element_text(size = 12, hjust = 0.5),
         axis.title.y = element_text(size = 12, hjust = 0.5),
         plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
         panel.background = element_rect(fill = "white"),
         legend.text = element_text(colour = "black", size = 10, hjust=0),
         legend.title = element_text(colour = "Black", size = 11,face = "bold", hjust = 0),
         legend.position = c(0.93, 0.19),
         panel.border = element_rect(colour = "gray",fill = NA, linewidth = 0.5),
         strip.text = element_text(size = 18, face = "italic", hjust = 0.5),
         axis.text = element_text(size = 12),
         legend.key.size = unit(0.5, "cm"),
         legend.spacing.y = unit(0.5, "cm"),
         legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgray"))

plot(volcanoPlot_NO_boxLabel)


                                      ### **************** ###
                                  ###   MAY BE SUPPLEMENTARY  ###
                                     ### **************** ###

############################################################################################
#####        6.          Differentially Expressed Proteins Only bar plot               #####
############################################################################################

# import data set 

DEPs_limma_output=read.csv("DEPs_limma_output.csv", fileEncoding = "latin1")

# visualize bar chart/plot
ggbarplot(DEPs_limma_output, x = "protein", y = "logFC",
          group = "Activity", fill = "Activity",           
          color = "white",            
          palette = c("blue", "red"),            
          sort.val = "desc",           
          sort.by.groups = T,     
          x.text.angle = 90,          
          ylab = "logFC",
          xlab = FALSE,
          legend.title = "Activity", 
          x.text.col = TRUE,                            
          ggtheme = theme_pubr())




############################################################################################
#####      7.        oxido-reductase activity proteins - GO Plot                       #####
############################################################################################


# import data set
oxidoreductase_GOterms <- read.csv("oxidoreductase_GOterms.csv", 
                                   row.names = 1, fileEncoding = "latin1")

# visualize circle plot
GOChord(oxidoreductase_GOterms,space = 0.04,gene.order = 'logFC',gene.space = 0.6,
        gene.size = 5,nlfc = 1,lfc.min = -3,lfc.max = 3,border.size = 0,process.label = 1,
        ribbon.col=c('#F8766D', '#7CAE00', '#00BFC4', '#C77CFF'),lfc.col=c('blue', 'red'))


############################################################################################
#####       8.       Plotting level-1 molecular function GO terms of DEPs               #####  
#####                              DEPs only                                           #####
############################################################################################

# import data set
GO_Terms_level_1 <- read.csv("GO_Terms_level_1.csv", encoding = "latin2") # data set containing only level1 GO terms

# reorder variables for facetting order
GO_Terms_level_1$Catagory <- factor(GO_Terms_level_1$Catagory , 
                                      levels = c("catalytic activity (GO:0003824)",
                                                 "binding (GO:0005488)",
                                                 "structural molecule activity (GO:0005198)",
                                                 "ATP-dependent activity (GO:0140657)"))

# represent GO term "yes", OR represent GO term "no"
# assign "" (empty) for proteins not represented in the specific specific GO term category ("no")
GO_Terms_level_1$protein[which(GO_Terms_level_1$GO_Term == "no")] <- ""

GO_Terms_level_1$protein <- iconv(GO_Terms_level_1$protein, from = "UTF-8", to = "ASCII//TRANSLIT")
GO_Terms_level_1$protein <- str_wrap(GO_Terms_level_1$protein, width = 20)

# Plot the results - bubble plot
buublePlot <- 
  ggplot(GO_Terms_level_1, aes(logFC, neg_log10_p_val)) +
  geom_point(aes(size = abundance_difference, 
                 color = GO_Term), 
             #shape = 21, 
             alpha = 0.5) +
  
  geom_text_repel(aes(label = protein, 
                      color = GO_Term), 
                  size = 2, 
                  alpha = 1,
                  #box.padding   = 0.15, 
                  #point.padding = 0.15, 
                  segment.color = 'darkblue',
                  segment.size = 0.1,
                  min.segment.length = 0.2,
                  lineheight = 0.8,
                  max.overlaps = 50) +
  
  scale_y_sqrt() +  # sqrt transformation of the y-axis
  scale_color_manual(values = c("darkgrey", "darkblue"), breaks = c("no", "yes"), labels = c("no", "yes")) +
  
  geom_vline(xintercept = 0, linetype = "solid", lwd = 0.1, color = "black") +  # Add vertical line at x = 0
  geom_hline(yintercept = c(-log10(0.05),-log10(0.01)), linetype = "solid", lwd = 0.1, color = "black") +  # Add vertical line at x = 0
  
  # Just for making it prettier!
  theme_light() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        text=element_text(size=12),
        plot.title = element_text(color = "black", 
                                  size = 12, 
                                  face = "bold", 
                                  hjust=1),
        strip.text = element_text(size = 9,
                                  face = "bold")) +
  
  guides(size = guide_legend(title = "Abundance difference  (%) between SAL and ALP",
                             color="black", 
                             hjust = 1, 
                             order = 1),
         fill = guide_legend(title="GO Term",
                             color="black", 
                             order = 2, 
                             size=10)) +
  
  facet_wrap(~Catagory, ncol = 2, labeller = label_wrap_gen(width=100))

options(repr.plot.width=14, repr.plot.height=14)
plot(buublePlot)


############################################################################################
#####         9.                    ALGAL AND FUNGAL genotypes                         #####
############################################################################################


##### import data set - algae_genotypes
algae_genotypes=read.csv("algae_genotypes.csv", fileEncoding = "latin1")
fungi_genotypes=read.csv("fungi_genotypes.csv", fileEncoding = "latin1")

## subset- algae
algae_genotypes_com <- algae_genotypes[3:ncol(algae_genotypes)]


#ordination by NMDS
nmds_algae_genotypes_com  = metaMDS(algae_genotypes_com, distance = "bray")

# data scores
data.scores_algae_nmds <- as.data.frame(scores(nmds_algae_genotypes_com)$sites)
data.scores_algae_nmds$SamplingSite = algae_genotypes$SamplingSite
data.scores_algae_nmds$Country = algae_genotypes$Country


# visualize NMDS plot of algal genotypes
algalGenotypes=ggplot(data.scores_algae_nmds, aes(x=NMDS1, y=NMDS2, color = SamplingSite))+ 
  geom_point(aes(shape = Country), cex=6) + 
  scale_color_manual(values=c("#228B22", "#00BFFF", "#00FFFF", 
                              "#FFD700", "#FF8C00", "#ADFF2F", 
                              "#B8860B", "#0000FF", "#A9A9A9", 
                              "#20B2AA", "#4682B4", "#F4A460")) + 
  scale_shape_manual(values=c(18, 15, 17, 19)) + theme_minimal()+
  theme (axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.title.x = element_text(size = 12, hjust = 0.5),
         axis.title.y = element_text(size = 12, hjust = 0.5),
         panel.background = element_rect(fill = "white"),
         legend.text = element_text(colour = "black", size = 10, hjust=0),
         legend.title = element_text(colour = "Black", size = 11,face = "bold", hjust = 0),
         legend.position = "right",
         panel.border = element_rect(colour = "gray",fill = NA, linewidth = 0.5),
         strip.text = element_text(size = 18,face = "italic", hjust = 0.5),
         axis.text = element_text(size = 12),
         legend.key.size = unit(0.5, "cm"),
         plot.title = element_text(hjust = 0.98, vjust = -10))+
  labs(title="Algal Genotypes")


plot(algalGenotypes)

## NMDS genotype FUNGI  
## subset- fungi
fungi_genotypes_com <- fungi_genotypes[3:ncol(fungi_genotypes)]


#ordination by NMDS
nmds_fungi_genotypes_com  = metaMDS(fungi_genotypes_com, distance = "bray")

# data scores
data.scores_fungi_nmds <- as.data.frame(scores(nmds_fungi_genotypes_com)$sites)
data.scores_fungi_nmds$SamplingSite = fungi_genotypes$SamplingSite
data.scores_fungi_nmds$Country = fungi_genotypes$Country


# visualize NMDS plot of fungal genotypes
fungalGenotypes=ggplot(data.scores_fungi_nmds, aes(x=NMDS1, y=NMDS2, color = SamplingSite))+ 
  geom_point(aes(shape = Country), cex=6) + 
  scale_color_manual(values=c("#228B22", "#00BFFF", "#00FFFF", 
                              "#FFD700", "#FF8C00", "#ADFF2F", 
                              "#B8860B", "#0000FF", "#A9A9A9", 
                              "#20B2AA", "#4682B4", "#F4A460")) + 
  scale_shape_manual(values=c(18, 15, 17, 19)) + theme_minimal()+
  theme (axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size = 11, face = "bold"),
         axis.title.x = element_text(size = 12, hjust = 0.5),
         axis.title.y = element_text(size = 12, hjust = 0.5),
         panel.background = element_rect(fill = "white"),
         legend.text = element_text(colour = "black", size = 10, hjust=0),
         legend.title = element_text(colour = "Black", size = 11,face = "bold", hjust = 0),
         legend.position = "right",
         panel.border = element_rect(colour = "gray",fill = NA, linewidth = 0.5),
         strip.text = element_text(size = 18,face = "italic", hjust = 0.5),
         axis.text = element_text(size = 12),
         legend.key.size = unit(0.5, "cm"),
         plot.title = element_text(hjust = 0.98, vjust = -10))+
  labs(title="Fungal Genotypes")

plot(fungalGenotypes)

# visualize combined plot
combined_genotype=ggarrange(fungalGenotypes, algalGenotypes, labels = c("A", "B"),
          common.legend = TRUE, legend = "right")

 plot(combined_genotype)
 