#!/usr/bin/Rscript
#This script is written for the microbiome analysis in the study "Effect of Nutritional Therapy on Microbiome in Canine enteropathy (¡®ENTiCE¡¯)".
#Written by Shuai Wang <wshuai@upenn.edu> in Dr. Beiting Lab <beiting@vet.upenn.edu> from University of Pennsylvania School of Veterinary Medicine.
#Date: Finished on 2018.10; Updated on 2019.3
#Note: The groups "Diet" and "Complex" in this script are equal to "DR" and "NDR" in the paper, respectively.
#Note: This script is written based on R version 3.4.2 (2017-09-28)
#Versions information:
sessionInfo()
# Versions of packages used in this script.
if (FALSE)
{
"R version 3.4.2 (2017-09-28)
attached base packages:
grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
WRS2_0.10-0         ggalluvial_0.6.0    btools_0.0.1        microbiomeSeq_0.1   corrplot_0.84       ggsci_2.8           ggord_1.0.0
phangorn_2.4.0      ape_5.0             devtools_1.13.5     ggpubr_0.1.6        magrittr_1.5        data.table_1.10.4-3 ggsignif_0.4.0
pheatmap_1.0.8      RColorBrewer_1.1-2  gplots_3.0.1        ggplot2_3.0.0       reshape2_1.4.3      scales_0.5.0        dplyr_0.7.4
vegan_2.4-6         lattice_0.20-35     permute_0.9-4       plyr_1.8.4          phyloseq_1.22.3

loaded via a namespace (and not attached):
Rcpp_1.0.0          mvtnorm_1.0-7       Biostrings_2.46.0   gtools_3.5.0        assertthat_0.2.0    digest_0.6.15       foreach_1.4.4
R6_2.2.2            stats4_3.4.2        pillar_1.2.1        zlibbioc_1.24.0     rlang_0.2.0         lazyeval_0.2.1      gdata_2.18.0
S4Vectors_0.16.0    Matrix_1.2-11       splines_3.4.2       stringr_1.3.0       igraph_1.2.1        munsell_0.4.3       compiler_3.4.2
pkgconfig_2.0.1     BiocGenerics_0.24.0 multtest_2.34.0     mgcv_1.8-20         biomformat_1.6.0    tibble_1.4.2        quadprog_1.5-5
IRanges_2.12.0      codetools_0.2-15    reshape_0.8.8       withr_2.1.2         MASS_7.3-47         bitops_1.0-6        mc2d_0.1-18
nlme_3.1-131        jsonlite_1.5        gtable_0.2.0        KernSmooth_2.23-15  stringi_1.1.7       XVector_0.18.0      bindrcpp_0.2
fastmatch_1.1-0     iterators_1.0.9     tools_3.4.2         ade4_1.7-10         Biobase_2.38.0      glue_1.2.0          parallel_3.4.2
survival_2.41-3     yaml_2.1.18         colorspace_1.3-2    rhdf5_2.22.0        cluster_2.0.6       caTools_1.17.1      memoise_1.1.0
bindr_0.1.1"
}


#*****************************************************Set working directory, load packages and define functions*********************************
setwd("working_path")
getwd()
#Loading packages
library("phyloseq")
library("plyr")
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library ("ggplot2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("ggsignif")
library("data.table")
library("ggpubr")
library(devtools)
library(phangorn)
library(ggord)
library("ggsci")
library("corrplot")
library("microbiomeSeq")
library("btools")
library(ggalluvial)
library("WRS2")
library("cowplot")
###################Define functions used in this script
#for mean, std  and error
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<- plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <-  plyr::rename(data_sum, c("mean" = varname))
 return(data_sum)
}
#trans
library(scales)
squish_trans <- function(from, to, factor,na.rm) {
  trans <- function(x,na.rm) {
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x,na.rm) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(trans_new("squished", trans, inv))
}


#Inputs
otutable <- import_biom(BIOMfilename = 'otu_table_mc2_w_tax_no_pynast_failures.json.biom',  treefilename = 'rep_set_pfiltered.midponit.fastree', parseFunction = parse_taxonomy_default)
mapping <- import_qiime_sample_data(mapfilename = 'ENTiCE_Mapping_file.txt')
physeq0 <- merge_phyloseq(otutable, mapping)
physeq0
rank_names(physeq0)

#Subsets of datasets for "DogFecal"
physeq_Dogfecal = subset_samples(physeq0, Description=="DogFecal")
physeq_Dogfecal
#Merge samples from different runs or sampleIDs
merged_physeq = merge_samples(physeq_Dogfecal,"SampleDescription", fun=mean)
#check
SD = merge_samples(sample_data(physeq_Dogfecal),"SampleDescription",fun=mean)
#Check here
identical(SD, sample_data(merged_physeq))
#sample_data(merged_physeq)$Owner <- levels(sample_data(physeq_Dogfecal)$Owner)[get_variable(merged_physeq,"Owner")]
sample_data(merged_physeq)$Timepoint2 <- levels(sample_data(physeq_Dogfecal)$Timepoint2)[get_variable(merged_physeq,
    "Timepoint2")]
sample_data(merged_physeq)$Remission <- levels(sample_data(physeq_Dogfecal)$Remission)[get_variable(merged_physeq,
    "Remission")]
sample_data(merged_physeq)$Breed_Strain <- levels(sample_data(physeq_Dogfecal)$Breed_Strain)[get_variable(merged_physeq,
    "Breed_Strain")]
sample_data(merged_physeq)$Sex<- levels(sample_data(physeq_Dogfecal)$Sex)[get_variable(merged_physeq,
    "Sex")]
sample_data(merged_physeq)$DogName<- levels(sample_data(physeq_Dogfecal)$DogName)[get_variable(merged_physeq,
    "DogName")]
sample_data(merged_physeq)$Description<- levels(sample_data(physeq_Dogfecal)$Description)[get_variable(merged_physeq,
    "Description")]
sample_data(merged_physeq)$SampleDescription<- levels(sample_data(physeq_Dogfecal)$SampleDescription)[get_variable(merged_physeq,
    "SampleDescription")]
sample_data(merged_physeq)$StudyID<- levels(sample_data(physeq_Dogfecal)$StudyID)[get_variable(merged_physeq,
    "StudyID")]

# check the samples
Day0Sample=subset(sample_data(merged_physeq), Timepoint2=="Day0")
sort(unique(Day0Sample$SampleDescription))
Day14Sample=subset(sample_data(merged_physeq), Timepoint2=="Day14")
sort(unique(Day14Sample$SampleDescription))
Day28Sample=subset(sample_data(merged_physeq), Timepoint2=="Day28")
sort(unique(Day28Sample$SampleDescription))
Day42Sample=subset(sample_data(merged_physeq), Timepoint2=="Day42")
sort(unique(Day42Sample$SampleDescription))
#
head(sample_data(merged_physeq))

#Taxonomic filtering
table(tax_table(merged_physeq)[, "Rank2"], exclude = NULL)
#The following ensures that features with ambiguous phylum annotation are also removed.
merged_physeq <- subset_taxa(merged_physeq, !is.na(Rank2) & !Rank2 %in% c("", "uncharacterized"))
table(tax_table(merged_physeq)[, "Rank2"], exclude = NULL)
# Remove Taxa without any reads
merged_physeq_f = prune_taxa(taxa_sums(merged_physeq) >=1, merged_physeq)
merged_physeq_f
#Remove singletons #Maybe sequencing errors.
merged_physeq_f  = filter_taxa(merged_physeq_f, function(x) sum(x>1)>=1, prune=TRUE)
#Verify taxonomic lineage names are correct and Change level names
colnames(tax_table(merged_physeq_f)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
tax_table(merged_physeq_f) <- gsub("D_\\d__","",tax_table(merged_physeq_f),perl=TRUE)
rank_names(merged_physeq_f)
#Check the reads for each sample
reads_each_sample = sample_sums(merged_physeq_f)
head(reads_each_sample)
#Rarefaction plot
#plot for rarefaction , run here only for the first time
#require("ranacapa")
#tiff(file="Rarefaction.tif", res = 600, width = 5500, height = 5800, compression = "lzw")
#p=ggrare(merged_physeq_f, step = 8, label = "SampleDescription", color = "Timepoint2", plot = TRUE, se = FALSE)
#p + facet_wrap(~Remission)
#dev.off()

#Get timepoints and groups for the analysis in this study
merged_physeq_f=subset_samples(merged_physeq_f, Remission%in%c("Diet","Complex","Healthy"))
merged_physeq_f1=subset_samples(merged_physeq_f, Timepoint2 %in% c("Day0", "Day14","Day28","Day42","Healthy"))
#Check the groups
sample_data(merged_physeq_f1)$Remission

#
merged_physeq_f1
#
rared_merged_physeq_f2= rarefy_even_depth (merged_physeq_f1,sample.size=10600, replace=FALSE, rngseed =750)



#Dogs with IBD at day 0
rared_merged_physeq_f2_overview = subset_samples(rared_merged_physeq_f2, Timepoint2 %in%c("Day0"))
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(rared_merged_physeq_f2_overview),
                 MARGIN = ifelse(taxa_are_rows(rared_merged_physeq_f2_overview), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_dataframe = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(rared_merged_physeq_f2_overview),
                      tax_table(rared_merged_physeq_f2_overview))
#Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf_dataframe, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Set up groups of "Diet" and "Complex"
rared_merged_physeq_f3  = subset_samples(rared_merged_physeq_f2, Remission %in% c("Diet","Complex"))

#*Plot disease scores
disease_scores = sample_data(rared_merged_physeq_f3)
disease_scores=subset(disease_scores, CCECAI !="NA")

#scores for all the dogs (complex +Diet)
disease_scores_day0 =subset(disease_scores, Timepoint2=="Day0")
disease_scores_day14 =subset(disease_scores, Timepoint2=="Day14")
mean(disease_scores_day0$CCECAI)
mean(disease_scores_day14$CCECAI)

#For "Diet' group
disease_scores_diet = subset(disease_scores, Remission=="Diet")
#Figure 1B
comparisons=list(c("Day0","Day14"),c("Day0","Day42"))
ggplot(disease_scores_diet, aes(Timepoint2,CCECAI,fill=Timepoint2))  + geom_violin(alpha = 0.9)+ylab("Abbrev. CCECAI of DR")+stat_compare_means(comparisons =comparisons,label = "p.signif",label.y.npc = "bottom")+
geom_boxplot(width=0.1, fill="white",outlier.shape = NA) +scale_fill_lancet()+scale_y_continuous(limits=c(0,9))+theme_bw()+
theme(axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 10))+theme(legend.position="none")
#Comparisons for "Diet' group
disease_scores_dietDay0 =subset(disease_scores_diet, Timepoint2=="Day0")
disease_scores_dietDay14 =subset(disease_scores_diet, Timepoint2=="Day14")
disease_scores_dietDay42 =subset(disease_scores_diet, Timepoint2=="Day42")
wilcox.test(disease_scores_dietDay0$CCECAI, disease_scores_dietDay14$CCECAI)
mean(disease_scores_dietDay0$CCECAI)
mean(disease_scores_dietDay14$CCECAI)
mean(disease_scores_dietDay42$CCECAI)

#Comparisons for "Complex" group
disease_scores_complex= subset(disease_scores, Remission=="Complex")
disease_scores_complexDay0 =subset(disease_scores_complex, Timepoint2=="Day0")
disease_scores_complexDay14 =subset(disease_scores_complex, Timepoint2=="Day14")
wilcox.test(disease_scores_complexDay0$CCECAI, disease_scores_complexDay14$CCECAI)
#Figure 1C
comparisons=list(c("Day28","Day42"),c("Day14","Day28"),c("Day0","Day14"),c("Day0","Day42"))
ggplot(disease_scores_complex, aes(Timepoint2,CCECAI,fill=Timepoint2)) + geom_violin(alpha = 0.9)+stat_compare_means(comparisons =comparisons,label = "p.signif",label.y.npc = "bottom")+
ylab("Abbrev. CCECAI of NDR")+geom_boxplot(width=0.1, fill="white",outlier.shape = NA) +
scale_fill_lancet()+scale_y_continuous(limits=c(0,21))+theme_bw()+
theme(axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 10))+theme(legend.position="none")

##################Compare the dogs with IBD and healthy dogs
rared_merged_physeq_f2

#Plot at phylum level
RelativeAb_physeq_bar = transform_sample_counts(rared_merged_physeq_f2, function(x) x / sum(x))
RelativeAb_physeq_bar_glom = tax_glom(RelativeAb_physeq_bar, taxrank="Phylum")
Relative_physeq_diet_bar_plot =data.table(psmelt(RelativeAb_physeq_bar_glom))
#Fig S2C
Relative_physeq_diet_bar_plot_pie = subset(Relative_physeq_diet_bar_plot, Timepoint2%in%c("Day0","Healthy"))
Relative_physeq_diet_bar_plot_pie$phylum=ifelse(Relative_physeq_diet_bar_plot_pie$Phylum%in%c("Fusobacteria", "Firmicutes","Actinobacteria", "Proteobacteria", "Bacteroidetes"), as.character(Relative_physeq_diet_bar_plot_pie$Phylum),"Others")
ggplot(Relative_physeq_diet_bar_plot_pie,aes(x="", y=Abundance,fill=phylum))+geom_bar(width = 1, stat = "identity",position = position_fill())+coord_polar("y")+theme_bw()+facet_wrap(~Timepoint2,nrow=2, strip.position="left")+
theme(text = element_text(size = 10),  legend.text=element_text(size=5), legend.box.spacing=unit(0,"line"), legend.margin=margin(0,0,0,0,'cm'), legend.key.height=unit(0.60,"line"),legend.key.width=unit(0.60,"line"),strip.text.x = element_blank(), panel.spacing.y = unit(-0.4,"lines"),panel.spacing.x = unit(-3, "lines"))+
theme(axis.text = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),strip.background= element_blank(),panel.border = element_blank())+
ylab("Relative Abundance")

#Compare Unifrac distances within groups #IBD vs. Healthy control.
RelativeAb_physeq_distance = subset_samples(RelativeAb_physeq_bar, Timepoint2 %in% c("Day0", "Healthy"))
DistunUif=phyloseq::distance(RelativeAb_physeq_distance, method="unifrac")
#day0 vs. helathy
DistunUif_martrix = as.matrix (DistunUif)
day0_distance=DistunUif_martrix[grep("_d0_",rownames(DistunUif_martrix), value=TRUE) ,grep("_d0_", colnames(DistunUif_martrix), value = TRUE)]
healthy_distance=DistunUif_martrix[grep("_dNA_", rownames(DistunUif_martrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_martrix), value = TRUE)]
day0_distance =subset(melt(day0_distance), Var1!=Var2)
healthy_distance =subset(melt(healthy_distance), Var1!=Var2)
library(rowr)
plot_dist=cbind.fill(as.vector(day0_distance$value),as.vector(healthy_distance$value), fill = NA)
colnames(plot_dist)= c("Day0","Healthy")
plot_dist2 = melt(plot_dist,variable.name="Timepoint")
my_comp= list(c("Day0","Healthy"))
#Fig S2D
ggplot(plot_dist2,aes(x=Timepoint, y =value, fill=Timepoint))+geom_boxplot(outlier.size =1)+stat_compare_means(comparisons=my_comp,label = "p.signif",label.y.npc = "top")+scale_fill_lancet()+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size=8), text = element_text(size = 10))+
theme(legend.position="none")+theme(axis.title.x = element_blank(),strip.background= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Unifrac distance within groups")+
scale_y_continuous(limits=c(NA,1.2))

#Phylogenetic distances
rared_merged_physeq_pd=subset_samples(rared_merged_physeq_f2, Timepoint2  %in% c("Day0", "Healthy"))
rared_merged_physeq_pd
pd <- estimate_pd(rared_merged_physeq_pd)
plot_pd= merge(pd, sample_data(rared_merged_physeq_pd), by=0, all=TRUE)
#Fig S2A
comparison = list(c("Day0","Healthy"))
ggplot(plot_pd, aes(x=Timepoint2,y=PD,fill=Timepoint2))+geom_boxplot(alpha=0.9,outlier.size =1) +scale_fill_lancet()+stat_compare_means(comparisons=my_comp,label = "p.signif")+
theme_bw()+scale_x_discrete("")+theme(axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size=8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 10))+
ylab("Phylogenetic Diversity")+guides(fill=FALSE)+scale_y_continuous(limits=c(NA,40))

#test
pd_h=subset(plot_pd,Timepoint2=="Healthy")
pd0=subset(plot_pd, Timepoint2=="Day0")
wilcox.test(pd0$PD,pd_h$PD)

#Shannon index
alpha_p2=plot_richness(rared_merged_physeq_pd,x="Timepoint2")
alphadt2 = data.table(alpha_p2$data)
# Subset to just the Shannon index
alpha_shannon <- alphadt2[(variable == "Shannon")]
# Order by Timepoint
alpha_shannon <- alpha_shannon[order(Timepoint2)][(is.na(se))]
scaleFUN <- function(x) sprintf("%.1f", x)
my_comparisons=list(c("Day0","Healthy"))
#Fig S2B
ggplot(alpha_shannon, aes(x=Timepoint2,y=value,group=Timepoint2,fill=Timepoint2))+geom_boxplot(alpha=0.9,outlier.size =1) +stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "bottom")+scale_fill_lancet()+theme_bw()+scale_x_discrete("")+
theme(axis.title.x=element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size=8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 10))+ylab("Shannon Index")+guides(fill=FALSE)+
scale_y_continuous(labels=scaleFUN, limits=c(NA,6.0))

#test
shannon_h =subset(alpha_shannon, Timepoint2=="Healthy")
shannon_IBD= subset(alpha_shannon, Timepoint2=="Day0")
wilcox.test(shannon_h$value,shannon_IBD$value)

#Filtering before ploting # maximum abundance >0.001 or in at least 10% samples
RelativeTernary_physeq = transform_sample_counts(rared_merged_physeq_f2, function(x) x / sum(x))
keepTaxa_Ternary= filter_taxa(RelativeTernary_physeq, function(x) max(x)>= 0.001, prune=FALSE)
Ternary_physeq_f= prune_taxa(keepTaxa_Ternary, rared_merged_physeq_f2)
Ternary_physeq_f2 =filter_taxa(Ternary_physeq_f, function(x) sum(x>0 )>=0.1*length(x),prune=TRUE)
Ternary_physeq_f3 = prune_taxa(taxa_sums(Ternary_physeq_f2) >=1, Ternary_physeq_f2)
#
physeq_ternary_f=subset_samples(Ternary_physeq_f3, Timepoint2%in%c("Day0","Healthy") &Remission%in%c("Diet","Complex","Healthy"))
melted_ternary_blom_OTUs  =psmelt(physeq_ternary_f)
melted_ternary_blom_OTUs$Remission =gsub("Diet", "DR", melted_ternary_blom_OTUs$Remission, perl=TRUE)
melted_ternary_blom_OTUs$Remission =gsub("Complex", "NDR", melted_ternary_blom_OTUs$Remission, perl=TRUE)
melted_phylum_top5 =subset (melted_ternary_blom_OTUs, Phylum %in% c("Fusobacteria", "Firmicutes","Actinobacteria", "Proteobacteria","Bacteroidetes"))
Phylum_abundance_ternary <- aggregate(Abundance~OTU+Remission+Phylum+Genus, melted_phylum_top5, FUN=mean)
Phylum_abundance_ternary_decast=dcast(Phylum_abundance_ternary, OTU+Phylum+Genus ~ Remission, value.var="Abundance")
Phylum_abundance_ternary_decast$Abundance =rowSums(Phylum_abundance_ternary_decast[,4:6])
#Figure 2A
library("ggtern")
ggtern(data = Phylum_abundance_ternary_decast, aes(x = DR, y = NDR, z =Healthy)) +geom_point(aes(fill =Phylum, size=log2(Abundance)),shape = 21)+scale_fill_npg()+theme_rgbw()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title = element_blank(), legend.position = c(0, 0.9), legend.justification = c(0, 1), legend.box.just = "left", legend.key=element_blank(), text=element_text(size=8),legend.text=element_text(size=8), legend.box.spacing=unit(0,"line"), legend.margin=margin(0,0,0,0,'cm'), legend.key.height=unit(0.3,"line"),legend.key.width=unit(0.3,"line"), panel.spacing.y = unit(-1,"lines"),panel.spacing.x = unit(-1, "lines"))+
guides(fill = guide_legend(override.aes = list(size=4)), size = "none")+scale_size_area(max_size=4)

#alpha diversity
#Get Phyloseq project for "Diet" group
rared_physeq_diet = subset_samples(rared_merged_physeq_f2, Remission %in%c("Diet"))

#phylogenetic distance (pd)
pd <- estimate_pd(rared_physeq_diet)
plot_pd= merge(pd, sample_data(rared_physeq_diet), by=0, all=TRUE)
#Figure S3A
scaleFUN <- function(x) sprintf("%.1f", x)
comparisons = list(c("Day0", "Day14"),c("Day0","Day28"),c("Day0","Day42"))
ggplot(plot_pd, aes(x=Timepoint2,y=PD,group=Timepoint2,fill=Timepoint2))+geom_boxplot(alpha=0.9,outlier.size =1)+scale_fill_lancet()+stat_compare_means(comparisons=comparisons, label = "p.signif")+
theme_bw()+scale_x_discrete("")+theme(axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
theme(text = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Phylogenetic Diversity")+guides(fill=FALSE)+
scale_y_continuous(labels=scaleFUN)
#test
pd0=subset(plot_pd, Timepoint2=="Day0")
pd14=subset(plot_pd, Timepoint2=="Day14")
wilcox.test(pd0$PD,pd14$PD)

#Shannon index for "diet' group
alpha_p2=plot_richness(rared_physeq_diet,x="Timepoint2",color="Remission")
alphadt2 = data.table(alpha_p2$data)
# Subset to just the Shannon index
alpha_shannon <- alphadt2[(variable == "Shannon")]
# Order by Timepoint
alpha_shannon <- alpha_shannon[order(Timepoint2)][(is.na(se))]
#For a paired comparison
alpha_shannon0=subset(alpha_shannon,Timepoint2=="Day0")
alpha_shannon14=subset(alpha_shannon,Timepoint2=="Day14")
wilcox.test(alpha_shannon0$value,alpha_shannon14$value)
#Try signed wilcox.test
#Figure S3B
comparisons = list(c("Day0", "Day14"),c("Day0","Day28"),c("Day0","Day42"))
ggplot(alpha_shannon, aes(x=Timepoint2,y=value,group=Timepoint2,fill=Timepoint2))+geom_boxplot(alpha=0.9,outlier.size =1) +scale_fill_lancet()+stat_compare_means(comparisons=comparisons, label = "p.signif")+
theme_bw()+scale_x_discrete("")+theme(axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
theme(text = element_text(size = 9),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Shannon Index")+guides(fill=FALSE)+
scale_y_continuous(labels=scaleFUN)

#Analyze based on dominant OTUs
RelativeAb_physeq_diet = transform_sample_counts(rared_physeq_diet, function(x) x / sum(x))
TopNOTUs <- names(sort(taxa_sums(RelativeAb_physeq_diet), TRUE)[1:40])
rared_physeq_diet_f   <- prune_taxa(TopNOTUs, rared_physeq_diet)
#Calculation of evenness
even_data =sample_data(rared_physeq_diet_f)
even_data$Shannon <- estimate_richness(rared_physeq_diet_f)[, "Shannon"]
even_data$Rich = specnumber(otu_table(rared_physeq_diet_f))
even_data$even = even_data$Shannon/log(even_data$Rich)
#Density of eveness for Top40
even_data_plot= subset(even_data, Timepoint2%in%c("Day0","Day14","Day28","Day42"))
#Figure 3A
ggplot(even_data_plot, aes(x=even, fill=Timepoint2)) +geom_density(alpha=0.7)+scale_x_continuous(limits=c(NA, 1))+
theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 9))+
scale_fill_lancet()

#qcomhd test for comparing the distributions of "day0" and "day14"
even_data_plot_test =subset(even_data_plot, Timepoint2%in%c("Day0","Day14"))
even_data_plot_melt= melt(even_data_plot_test, id.vars=c("Timepoint2","even"))
#significant; P-value <0.05 around the peaks
qcomhd(even~Timepoint2, data=even_data_plot_melt, q = seq(0, 1, by=.1))

##The phylogenetic distance to healthy dogs
rared_merged_physeq_dis=subset_samples(rared_merged_physeq_f2, Remission  %in% c("Diet", "Healthy"))
Relative_physeq_dis  = transform_sample_counts(rared_merged_physeq_dis, function(x) x / sum(x))
DistunUif=phyloseq::UniFrac(Relative_physeq_dis,weighted=FALSE)
DistunUif_matrix= as.matrix(DistunUif)
#days
day0_distance=DistunUif_matrix[grep("_d0_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day14_distance=DistunUif_matrix[grep("_d14_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day28_distance=DistunUif_matrix[grep("_d28_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day42_distance=DistunUif_matrix[grep("_d42_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
library(rowr)
plot_dist=cbind.fill(as.vector(day0_distance),as.vector(day14_distance),as.vector(day28_distance),as.vector(day42_distance), fill = NA)
colnames(plot_dist)= c("Day0","Day14","Day28","Day42")
plot_dist2 = melt(plot_dist,variable.name="Timepoint")
#Figure 3C
my_comparisons = list (c("Day0","Day14"),c("Day0","Day28"),c("Day0","Day42"))
ggplot(plot_dist2, aes(x=Timepoint, y =value,fill=Timepoint))+geom_boxplot(outlier.size=1)+stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "bottom")+scale_fill_lancet() +
ylab("Unifrac distance")+theme_bw()+theme(legend.position="none", strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 9))+
scale_y_continuous(limits=c(NA,1.1))
wilcox.test(plot_dist$Day0, plot_dist$Day42)


#Fig3D for alluvium plot for "Diet"
RelativeAb_physeq_diet = transform_sample_counts(rared_physeq_diet, function(x) x / sum(x))
rared_physeq_diet_tax_glom_phylum = tax_glom(RelativeAb_physeq_diet, taxrank="Phylum")
rared_physeq_diet_melted_phylum  =psmelt(rared_physeq_diet_tax_glom_phylum)
rared_physeq_diet_melted_phylum_top5 =subset (rared_physeq_diet_melted_phylum, Phylum %in% c("Fusobacteria", "Firmicutes","Actinobacteria", "Proteobacteria","Bacteroidetes"))
Phylum_abundance <- aggregate(Abundance~Timepoint2+Phylum, rared_physeq_diet_melted_phylum_top5, FUN=mean)
ggplot(data=Phylum_abundance , aes(x=Timepoint2, weight=Abundance, alluvium=Phylum))+geom_alluvium(aes(fill=Phylum), alpha=0.75, decreasing =FALSE)+theme_minimal()+xlab("Timepoint in the trial")+
theme(legend.text=element_text(size=6), legend.position="top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(12,0,-10,0), legend.key.width=unit(0.4,'line'),legend.key.height=unit(0.4,'line'),legend.title = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 9))



#Alpha diversity
#Get phyloseq project for "Complex"
rared_merged_physeq_complex  = subset_samples(rared_merged_physeq_f2, Remission =="Complex")
all_even_data =sample_data(rared_merged_physeq_complex)
all_even_data$Shannon <- estimate_richness(rared_merged_physeq_complex )[, "Shannon"]
all_even_data$Rich = specnumber(otu_table(rared_merged_physeq_complex))
all_even_data$even = all_even_data$Shannon/log(all_even_data$Rich)
#Figure 4A
even_data_plot_all= subset(all_even_data, Timepoint2%in%c("Day0","Day14","Day28", "Day42"))
ggplot(even_data_plot_all, aes(x=even, fill=Timepoint2)) +geom_density(alpha=0.7)+scale_x_continuous(limits=c(NA, 1))+
theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 9))+
scale_fill_lancet()
#test
even_data_plot_test =subset(even_data_plot_all, Timepoint2%in%c("Day0","Day14"))
even_data_plot_all_melt= melt(even_data_plot_test, id.vars=c("Timepoint2","even"))
qcomhd(even~Timepoint2, data=even_data_plot_all_melt, q = seq(0, 1, by=.1))


##the distance between healthy dogs and "Complex"
rared_merged_physeq_dis=subset_samples(rared_merged_physeq_f2, Remission  %in% c("Complex", "Healthy"))
Relative_physeq_dis  = transform_sample_counts(rared_merged_physeq_dis, function(x) x / sum(x))
DistunUif=phyloseq::UniFrac(Relative_physeq_dis,weighted=FALSE)
DistunUif_matrix= as.matrix(DistunUif)
#day0
day0_distance=DistunUif_matrix[grep("_d0_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day14_distance=DistunUif_matrix[grep("_d14_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day28_distance=DistunUif_matrix[grep("_d28_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
day42_distance=DistunUif_matrix[grep("_d42_",rownames(DistunUif_matrix), value=TRUE) ,grep("_dNA_", colnames(DistunUif_matrix), value = TRUE)]
library(rowr)
plot_dist=cbind.fill(as.vector(day0_distance),as.vector(day14_distance),as.vector(day28_distance),as.vector(day42_distance), fill = NA)
colnames(plot_dist)= c("Day0","Day14","Day28","Day42")
plot_dist2 = melt(plot_dist,variable.name="Timepoint")
#Figure 4B
my_comparisons = list (c("Day0","Day14"),c("Day0","Day28"),c("Day14","Day28"), c("Day28","Day42"))
ggplot(plot_dist2, aes(x=Timepoint, y =value,fill=Timepoint))+geom_boxplot(outlier.size=1)+stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "bottom")+
ylab("Unifrac distance")+scale_fill_lancet() +theme_bw()+theme(axis.title.x=element_blank(), legend.position="none",axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 9))+
scale_y_continuous(limits=c(NA,1.2) )
wilcox.test(plot_dist$Day0, plot_dist$Day42)


#Figure 4C #alluvium
Re_rared_merged_physeq_f2 = transform_sample_counts(rared_merged_physeq_f2,function(x) x / sum(x))
Re_rared_merged_physeq_f2_complex = subset_samples(Re_rared_merged_physeq_f2, Remission =="Complex")
rared_physeq_complex_tax_glom_phylum = tax_glom(Re_rared_merged_physeq_f2_complex, taxrank="Phylum")
rared_physeq_complex_melted_phylum  =psmelt(rared_physeq_complex_tax_glom_phylum)
rared_physeq_complex_melted_phylum_top5 =subset (rared_physeq_complex_melted_phylum, Phylum %in% c("Fusobacteria", "Firmicutes","Actinobacteria", "Proteobacteria","Bacteroidetes"))
Phylum_abundance <- aggregate(Abundance~Timepoint2+Phylum, rared_physeq_complex_melted_phylum_top5, FUN=mean)
ggplot(data=Phylum_abundance , aes(x=Timepoint2, weight=Abundance, alluvium=Phylum))+geom_alluvium(aes(fill=Phylum), alpha=0.75, decreasing =FALSE)+theme_minimal()+ylab("Relative Abundance")+
theme(axis.title.x=element_blank(), legend.text=element_text(size=6), legend.position="top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(12,0,-10,0), legend.key.width=unit(0.4,'line'),legend.key.height=unit(0.4,'line'),legend.title = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 9))

#Beta diversity
#comparison between healthy dogs and dogs with IBD
Relative_physeq_diet_beta  = transform_sample_counts(rared_physeq_diet, function(x) x / sum(x))
Relative_physeq_diet = subset_samples(Relative_physeq_diet_beta, Timepoint2 %in% c("Day0", "Day14","Day42"))
DistunUif=phyloseq::distance(Relative_physeq_diet ,method="wunifrac")
ordunUif = ordinate(Relative_physeq_diet , method = "PCoA", distance = DistunUif)
plot_ordination(Relative_physeq_diet , ordunUif,color="DogName") + ggtitle("PCoA: weighted Unifrac")+ geom_path()+
geom_point(aes(shape=factor(Timepoint2)),size=5)+scale_shape_manual(name="TimePoint",values = c('Day0' = 17, 'Day14' = 16,'Day28'=15,'Day42'=14))+
theme_bw() + theme(text = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color=FALSE)

#Figure 3B
DistunUif=phyloseq::distance(Relative_physeq_diet ,method="unifrac")
ordunUif = ordinate(Relative_physeq_diet , method = "PCoA", distance = DistunUif)
plot_ordination(Relative_physeq_diet, ordunUif, color="DogName") + geom_path(lineend=NULL,linejoin=NULL)+scale_fill_lancet() +
geom_point(aes(shape=factor(Timepoint2)),size=2)+scale_shape_manual(name="TimePoint",values = c('Day0' = 17, 'Day14' = 16, 'Day42'=15))+
theme_bw()+guides(color=FALSE)+
theme(legend.title = element_blank(),legend.key.height=unit(0.49,"line") ,legend.key.width=unit(0.5,"line"), legend.position=c(0.8,0.15),legend.text=element_text(size=5), legend.box.spacing=unit(0,"line"), legend.margin=margin(0,0,0,0,'cm'))+
theme(text = element_text(size = 9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+guides(fill=FALSE)



#Set up a new name of 'merged_physeq_f1'
merged_physeq_differ = merged_physeq_f1

# Library packages #This may mask some functions from other packages
library("DESeq2")
packageVersion("DESeq2")

#Set up comparison between dogs with IBD and healthy dogs
merged_physeq_f_deseq= prune_samples(names(which(sample_sums(merged_physeq_differ) >= 10000)), merged_physeq_differ)
#filtering # maximum abundance >0.001 or inat least 10% samples
RelativeDeSeq_physeq = transform_sample_counts(merged_physeq_f_deseq, function(x) x / sum(x))
keepTaxa_Deseq= filter_taxa(RelativeDeSeq_physeq, function(x) max(x)>= 0.001, prune=FALSE)
Deseq_physeq_f= prune_taxa(keepTaxa_Deseq, merged_physeq_f_deseq)
Deseq_physeq_f2 =filter_taxa(Deseq_physeq_f, function(x) sum(x>0 )>=0.1*length(x),prune=TRUE)
Deseq_physeq_f3 = prune_taxa(taxa_sums(Deseq_physeq_f2) >=1, Deseq_physeq_f2)
contrast_day0vsHealthy = subset_samples(Deseq_physeq_f3, Timepoint2 %in% c("Day0", "Healthy"))
#OTU level
deseqstudy = phyloseq_to_deseq2(contrast_day0vsHealthy, ~Timepoint2)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res = results(deseqstudy, cooksCutoff = FALSE,pAdjustMethod = "BH",independentFiltering=FALSE)
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(contrast_day0vsHealthy)[rownames(sigtab), ], "matrix"))
head(sigtab)
#Supplementary Table S1
write.table(sigtab, "DESeq_HealthyVS.Day0_OTUs_sig.xls", sep="\t")

#For 'Diet' Group
#Choose group to be compared
merged_physeq_f_deseq= prune_samples(names(which(sample_sums(merged_physeq_differ) >= 10000)), merged_physeq_differ)
#Get Diet.
Deseq_physeq_diet= subset_samples(merged_physeq_f_deseq, Remission=="Diet")
#filtering # maximum abundance >0.001 or inat least 10% samples
RelativeDeSeq_physeq_diet = transform_sample_counts(Deseq_physeq_diet, function(x) x / sum(x))
keepTaxa_Deseq= filter_taxa(RelativeDeSeq_physeq_diet, function(x) max(x)>= 0.001, prune=FALSE)
Deseq_physeq_diet_f= prune_taxa(keepTaxa_Deseq, Deseq_physeq_diet)
Deseq_physeq_diet_f2 =filter_taxa(Deseq_physeq_diet_f, function(x) sum(x>0 )>=0.1*length(x),prune=TRUE)
Deseq_physeq_diet_f2 = prune_taxa(taxa_sums(Deseq_physeq_diet_f2) >=1, Deseq_physeq_diet_f2)
contrast_day0vs14 = subset_samples(Deseq_physeq_diet_f2, Timepoint2 %in% c("Day0", "Day14"))
#<strong>For a species level</strong>
deseqstudy = phyloseq_to_deseq2(contrast_day0vs14, ~Timepoint2)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res_species = results(deseqstudy, cooksCutoff = FALSE,pAdjustMethod = "BH",independentFiltering=FALSE)
alpha = 0.05
sigtab = res_species[which(res_species$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(contrast_day0vs14)[rownames(sigtab), ], "matrix"))
head(sigtab)
#Supplementary Table S4
write.table(sigtab, "DESeq_Diet_Day14.VS.Day0.xls", sep="\t")

#Valcano plot
valcano_species = res_species
valcano_species$LogPvalue = -log10(valcano_species$pvalue)
valcano_species_diff = cbind(as(valcano_species, "data.frame"), as(tax_table(contrast_day0vs14)[rownames(valcano_species), ], "matrix"))
valcano_species_diff <-valcano_species_diff %>% plyr::mutate(Differential=ifelse(LogPvalue  >1.301 & log2FoldChange >=1, "Up", ifelse(LogPvalue  >1.301 & log2FoldChange <=-1,"Down", "NoChange")))
library(ggrepel)
ggplot(valcano_species_diff , aes(x = log2FoldChange, y =LogPvalue,color=Differential))+geom_point()+scale_color_manual(values=c("#E64B35FF","grey70","#4DBBD5FF"))+ scale_shape_manual(values=c(0, 1, 16,10))+
geom_hline(yintercept = 1.301, colour="#990000", linetype="dashed") + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed")+
ylab("-log10Pvalue")+theme_bw()+theme(legend.position="none", panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
theme(text = element_text(size = 9),plot.title = element_text(size = 8, hjust = 0.5,vjust=0))+ggtitle("Days 0 versus 14")+scale_y_continuous(labels=scaleFUN)+
geom_text_repel(aes(label=ifelse(Genus%in%c("Clostridium sensu stricto 1", "Escherichia-Shigella","Blautia","Megasphaera") & Differential !="NoChange",as.character(Genus),'')),hjust=0,vjust=0,size=2.2)


#For a genus level #"Diet" group
physeqGenus = tax_glom(contrast_day0vs14, "Genus")
deseqstudy = phyloseq_to_deseq2(physeqGenus, ~Timepoint2)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res = results(deseqstudy, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeqGenus)[rownames(sigtab), ], "matrix"))
head(sigtab)
#Supplementary Table S3
write.table(sigtab, "DESeq_Diet_genus_Days14_vs_0.xls", sep="\t")


#*****For comlex group
#Choose group to be compared
merged_physeq_f_deseq= prune_samples(names(which(sample_sums(merged_physeq_differ) >= 10000)), merged_physeq_differ)
Deseq_physeq_cox= subset_samples(merged_physeq_f_deseq, Remission=="Complex")
#filtering # maximum abundance >0.001 or inat least 10% samples
RelativeDeSeq_physeq_cox = transform_sample_counts(Deseq_physeq_cox, function(x) x / sum(x))
keepTaxa_Deseq= filter_taxa(RelativeDeSeq_physeq_cox, function(x) max(x)>= 0.001, prune=FALSE)
Deseq_physeq_cox_f= prune_taxa(keepTaxa_Deseq, Deseq_physeq_cox)
Deseq_physeq_cox_f2 =filter_taxa(Deseq_physeq_cox_f, function(x) sum(x>0)>=0.1*length(x),prune=TRUE)
Deseq_physeq_cox_f2 = prune_taxa(taxa_sums(Deseq_physeq_cox_f2) >=1, Deseq_physeq_cox_f2)
contrast_day0vs14 = subset_samples(Deseq_physeq_cox_f2, Timepoint2 %in% c("Day0", "Day14"))
#For a species level
deseqstudy = phyloseq_to_deseq2(contrast_day0vs14, ~Timepoint2)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res_species = results(deseqstudy, cooksCutoff = FALSE,pAdjustMethod = "BH",independentFiltering=FALSE)
alpha = 0.05
sigtab = res_species[which(res_species$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(contrast_day0vs14)[rownames(sigtab), ], "matrix"))
head(sigtab)
#Table S5
write.table(sigtab, "DESeq_complex_Days14_vs_0.xls", sep="\t")

#For a genus level #The ID is randomly selected as a representive, so not the exact one.
physeqGenus = tax_glom(contrast_day0vs14, "Genus")
deseqstudy = phyloseq_to_deseq2(physeqGenus, ~Timepoint2)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res = results(deseqstudy, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeqGenus)[rownames(sigtab), ], "matrix"))
head(sigtab)
#Table S4
write.table(sigtab, "DESeq_complex_genus_Days14_vs_0.xls", sep="\t")


#Comparison of 'complex' and diet groups at day 0
# Library packages
library("DESeq2")
#Choose group to be compared
merged_physeq_f_deseq= prune_samples(names(which(sample_sums(merged_physeq_differ) >= 10000)), merged_physeq_differ)
Deseq_physeq_coxVSdiet= subset_samples(merged_physeq_f_deseq, Remission %in% c("Complex","Diet"))
#filtering # maximum abundance >0.001 or inat least 10% samples
RelativeDeSeq_physeq_dietVScox = transform_sample_counts(Deseq_physeq_coxVSdiet, function(x) x / sum(x))
keepTaxa_Deseq= filter_taxa(RelativeDeSeq_physeq_dietVScox, function(x) max(x)>= 0.001, prune=FALSE)
Deseq_physeq_coxVSdiet_f= prune_taxa(keepTaxa_Deseq, Deseq_physeq_coxVSdiet)
Deseq_physeq_coxVSdiet_f2 =filter_taxa(Deseq_physeq_coxVSdiet_f, function(x) sum(x>0 )>=0.1*length(x),prune=TRUE)
Deseq_physeq_coxVSdiet_f2 = prune_taxa(taxa_sums(Deseq_physeq_coxVSdiet_f2) >=1, Deseq_physeq_coxVSdiet_f2)
contrast_coxvsdiet = subset_samples(Deseq_physeq_coxVSdiet_f2, Timepoint2 %in% c("Day0"))
# export a tree file for Graphlan
tree1 = phy_tree(contrast_coxvsdiet)
ape::write.tree(tree1, "contrast_coxvsdiet.tree")
#For a species level
deseqstudy = phyloseq_to_deseq2(contrast_coxvsdiet, ~Remission)
#tolerate 0
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqstudy), 1, gm_mean)
deseqstudy = estimateSizeFactors(deseqstudy, geoMeans = geoMeans)
deseqstudy = DESeq(deseqstudy, test="Wald", fitType="parametric")
res_species = results(deseqstudy, cooksCutoff = FALSE,pAdjustMethod = "BH",independentFiltering=FALSE)
alpha = 0.05
sigtab = res_species[which(res_species$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(contrast_coxvsdiet)[rownames(sigtab), ], "matrix"))
#Table S2
write.table(sigtab, "DESeq_Complex_vs_diet_dya0.xls", sep="\t")

#Figure 2F
head(sigtab)
library("ggrepel")
plot_fig = data.frame(sigtab, OTUs = rownames(sigtab))
plot_fig=plot_fig[order(plot_fig$Phylum),]
plot_fig$OTUs = factor(plot_fig$OTUs, levels=plot_fig$OTUs)
ggplot(data = plot_fig,
  # group lets `ggplot` know we want different errorbars/bars for each day
       aes(x = OTUs, y = `log2FoldChange`, group = Phylum)
    ) +
    geom_bar(
        stat = "identity",
        aes(fill = Phylum),
        position = position_dodge(width = 0.9)
    ) +
    ylab("Fold Change (log2)")+scale_fill_npg() +theme_bw()+
    theme(axis.text.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 8))+xlab("OTUs with differential abundance")+
    theme(legend.margin=margin(t = 0, unit='cm'), legend.key.height=unit(0.8,"line"),legend.key.width=unit(0.8,"line"))+geom_text_repel(aes(label=ifelse(Genus%in%c("Clostridium sensu stricto 1", "Escherichia-Shigella") ,as.character(Genus),'')),hjust=0,vjust=0,size=2.2)


#Comparison of genus with diferential abundance
comparison_of_diet_complex_response =read.table(file="comparison_of_diet_complex_response.txt", header=T,sep="\t")
comparison_of_diet_complex_response$Genus=factor(comparison_of_diet_complex_response$Genus,levels =rev(levels(comparison_of_diet_complex_response$Genus)))
#Figure 4D
ggplot(comparison_of_diet_complex_response, aes(x=Group, y=Genus)) + geom_point(aes(size=Abs_Log2foldchange,color=Difference))+scale_size_continuous(range = c(1,4))+scale_color_lancet()+
theme_bw()+theme(legend.margin=margin(t = 0, unit='cm'), legend.key.height=unit(0.6,"line"),legend.key.width=unit(0.6,"line"), axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
theme(legend.text=element_text(size=6), legend.title=element_text(size=7), legend.box.spacing=unit(0.1,"line"), legend.margin=margin(0,0,0,0,'cm'), legend.key.height=unit(0.60,"line"),legend.key.width=unit(0.60,"line"))+
theme(axis.text.y=element_text(face="italic"), axis.title.x=element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 9))+scale_x_discrete(labels = c("DR","NDR"))



#Get phyloseq for Associations
rared_physeq_all_AC = subset_samples(rared_merged_physeq_f2, Remission %in%c("Diet","Complex"))
#Mix the lables of Genus and species for plot # Define the ranks you want to include
myranks = c("Genus", "Species")
mylabels = apply(tax_table(rared_physeq_all_AC)[, myranks], 1, paste, sep="", collapse="_")
# Add concatenated labels as a new rank after strain
tax_table(rared_physeq_all_AC) <- cbind(tax_table(rared_physeq_all_AC), catglab=mylabels)
#To test the association for each species, we first filtered low abundance species and
#confined our analysis to XX species that accounted for at least 0.01% of microbial composition and were present in more than 10% participants
RelativeCorr_physeq_all_AC = transform_sample_counts(rared_physeq_all_AC, function(x) x / sum(x))
keepTaxa_Corr= filter_taxa(RelativeCorr_physeq_all_AC, function(x) max(x)>= 0.001, prune=FALSE)
Corr_physeq_all_f= prune_taxa(keepTaxa_Corr, rared_physeq_all_AC)
physeq_all_corr =filter_taxa(Corr_physeq_all_f, function(x) sum(x>0 )>=0.1*length(x),prune=TRUE)
physeq_all_corr = prune_taxa(taxa_sums(physeq_all_corr) >=1, physeq_all_corr)
#Check
physeq_all_corr
identical(rownames((otu_table(physeq_all_corr))),rownames(sample_data(physeq_all_corr)))
#Start analysis
require (microbiomeSeq)
## Log-ransformed of rarefied OTU table
LogTransformed_physeq_all = transform_sample_counts(physeq_all_corr, function(x) log(x+1))
#Remove the phylum with otu number nearly 0 and transform the non-meanningful variables into 'NA'
physeq_all_correlation <- subset_taxa(LogTransformed_physeq_all , !(Phylum %in% c("Planctomycetes","Verrucomicrobia","Euryarchaeota","Spirochaetae")))
#Remove the Viriables that contain too many 'NA' or useless
sample_data(physeq_all_correlation)$X.SampleID <- NA
sample_data(physeq_all_correlation)$Heptanoic_Acid <-NA
sample_data(physeq_all_correlation)$Caproic_Acid <-NA
sample_data(physeq_all_correlation)=sample_data(physeq_all_correlation)[, colSums(sample_data(physeq_all_correlation) != 0, na.rm = TRUE) > 5]
#Genus
    physeq_all_cor_genus <- taxa_level(physeq_all_correlation, "Genus")
env.taxa.cor_genus<- taxa.env.correlation(physeq_all_cor_genus, grouping_column = "Description", method = "spearman",
    pvalue.threshold = 0.05, padjust.method = "fdr", adjustment = 5, num.taxa = 10000)
    write.table(env.taxa.cor_genus, "Diet_complex_corr_Genus.xls", sep="\t")
    #Species
    #Take care of the names here. Results will be wrong if you use a "Species" as a selected level.
    physeq_all_cor_species <- taxa_level(physeq_all_correlation, "catglab")
env.taxa.cor_species <- taxa.env.correlation(physeq_all_cor_species, grouping_column = "Description", method = "spearman",
    pvalue.threshold = 0.05, padjust.method = "fdr", adjustment = 5, num.taxa = 10000, select.variables = NULL)
#Table S8
    write.table(env.taxa.cor_species, "Diet_complex_corr_species.xls", sep="\t")

#Plot correlation for Ecoli
slected_OTU_Ecoli = subset(otu_table(LogTransformed_physeq_all), select= c("FJ950694.1.1472"))
plot_OTU_cor_Ecoli = merge(slected_OTU_Ecoli , sample_data(LogTransformed_physeq_all), by=0, all=TRUE)
#Figure 2C
melted_plot_otu_cor_Ecoli<- melt(plot_OTU_cor_Ecoli, id.vars=c("CCECAI","DogName","Timepoint2"),measure.vars="FJ950694.1.1472",variable.name="OTU")
ggplot(melted_plot_otu_cor_Ecoli, aes(CCECAI,value)) +  geom_jitter(alpha=1,size=0.8) +   stat_smooth(method=lm) +ylab("Log(Abundance)\nE. coli")+theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 8))+xlab("Disease scores")+
theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(strip.background = element_blank(),plot.title = element_text(size = 8,vjust=0,hjust=0.5))
#test
cor.test(plot_OTU_cor_Ecoli$CCECAI, plot_OTU_cor_Ecoli$FJ950694.1.1472,method="spearman")


#Plot for C. perfringens
slected_OTU_clos = subset(otu_table(LogTransformed_physeq_all), select= c("New.ReferenceOTU52"))
plot_OTU_cor_clos = merge(slected_OTU_clos, sample_data(LogTransformed_physeq_all), by=0, all=TRUE)
#Figure 2E
melted_plot_otu_cor_clos<- melt(plot_OTU_cor_clos, id.vars=c("CCECAI","DogName","Timepoint2"),measure.vars="New.ReferenceOTU52",variable.name="OTU")
ggplot(melted_plot_otu_cor_clos, aes(CCECAI,value)) +  geom_jitter(alpha=1,size=0.8) +   stat_smooth(method=lm) +ylab("Log(Abundance)\nClostridium sp.")+theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 8))+xlab("Disease scores")+
theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(strip.background = element_blank(),plot.title = element_text(size = 8,vjust=0,hjust=0.5))
#test
cor.test(plot_OTU_cor_clos$CCECAI,plot_OTU_cor_clos$New.ReferenceOTU52,method="spearman")



merged_physeq_differ = merged_physeq_f1
#********Phyloseq used# Relative abundance
RA_merged_physeq_differ = transform_sample_counts(merged_physeq_differ, function(x) x/sum(x))
#*************************

#Relative abundance #Healthy VS. disease
RA_merged_physeq_differ_f=subset_samples(RA_merged_physeq_differ, Timepoint2%in%c("Day0","Healthy")&Remission%in%c("Diet","Complex","Healthy"))

##For E .coli
physeq_Ecoli <- subset_taxa(RA_merged_physeq_differ_f, rownames(tax_table(RA_merged_physeq_differ_f))=="FJ950694.1.1472")
sub_beta_RA_data <- data.table(psmelt(physeq_Ecoli))
sub_beta_RA_data$Timepoint2 = factor(sub_beta_RA_data$Timepoint2, levels =c("Day0", "Healthy"))
scaleFUN <- function(x) sprintf("%.1f", x)
my_comp =list(c("Day0","Healthy"))
ggplot(sub_beta_RA_data, aes(x=Timepoint2,y=Abundance,fill=Timepoint2))+geom_boxplot(outlier.size =0.6)+stat_compare_means(comparisons=my_comp,label = "p.signif",label.y.npc = "top")+scale_fill_manual(values=c("#F39B7FFF","#F39B7FFF"))+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size=8), text = element_text(size = 8))+
theme(legend.position="none")+theme(axis.title.x = element_blank(),strip.background= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance\nof E. coli")+
scale_y_continuous(limits=c(0,1.4), breaks=seq(0,1,by=0.2),trans=squish_trans(0,0.2,0.1))



#For C. perfringens
physeq_clostri_com <- subset_taxa(RA_merged_physeq_differ, rownames(tax_table(RA_merged_physeq_differ))=="New.ReferenceOTU52")
physeq_clostri_com_f=subset_samples(physeq_clostri_com, Timepoint2%in%c("Day0","Healthy")&Remission%in%c("Diet","Complex","Healthy"))
sub_beta_RA_com_f <- data.table(psmelt(physeq_clostri_com_f))
#Figure 2D
ggplot(sub_beta_RA_com_f, aes(x=Timepoint2, y=Abundance,fill=Timepoint2))+ geom_boxplot(outlier.size =0.6) +ylab("Relative Abundance\nof Clostridium sp.")+stat_compare_means(label="p.signif", comparisons=list(c("Day0","Healthy")))+
scale_fill_manual(values=c("#00A087FF","#00A087FF"))+theme_bw()+theme(legend.position="none", axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 8))+
scale_y_continuous(limits=c(0,1.4), breaks=seq(0,1,by=0.2),trans=squish_trans(0,0.2,0.1))

#All time points for "Diet" and healthy Groups
RA_merged_physeq_differ_f=subset_samples(RA_merged_physeq_differ, Remission %in% c("Diet","Healthy"))
physeq_Ecoli <- subset_taxa(RA_merged_physeq_differ_f, rownames(tax_table(RA_merged_physeq_differ_f))=="FJ950694.1.1472")
sub_beta_RA_data <- data.table(psmelt(physeq_Ecoli))
sub_beta_RA_data$Timepoint2 = factor(sub_beta_RA_data$Timepoint2, levels =c("Day0", "Day14","Day28","Day42","Healthy"))
#Figure 3F# paired comparison
ggplot(sub_beta_RA_data, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =1) +
geom_signif(annotations=c("*","ns"), y_position=c(11.0,11.5), xmin=c(1,4), xmax=c(4,5),size=0.3,textsize=3)+
ylab("Relative Abundance\nof E.coli")+scale_fill_lancet()+theme_bw()+
theme(axis.title.x=element_blank(), legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,3), breaks=c(0, 0.2,0.6,1.0), trans=squish_trans(0,0.1,0.01))


#C. perfringens in "Diet" and healthy
RA_merged_physeq_differ_f=subset_samples(RA_merged_physeq_differ, Remission %in% c("Diet","Healthy"))
physeq_clostri <- subset_taxa(RA_merged_physeq_differ_f, rownames(tax_table(RA_merged_physeq_differ_f))=="New.ReferenceOTU52")
sub_beta_RA_data <- data.table(psmelt(physeq_clostri))
sub_beta_RA_data$Timepoint2 = factor(sub_beta_RA_data$Timepoint2, levels =c("Day0", "Day14","Day28","Day42","Healthy"))
#Figure 3G
ggplot(sub_beta_RA_data, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =1) +ylab("Relative Abundance\nof Clostridium sp.")+stat_compare_means(label="p.signif", comparisons=list(c("Day0","Day42"),c("Day42","Healthy")))+
scale_fill_lancet()+theme_bw()+theme(axis.title.x=element_blank(),legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,3), breaks=c(0, 0.2,0.6,1.0),trans=squish_trans(0,0.1,0.01))


#C. hiranonis in Diet
RA_merged_physeq_differ_f=subset_samples(RA_merged_physeq_differ, Remission %in% c("Diet"))
physeq_hiranonis <- subset_taxa(RA_merged_physeq_differ_f, rownames(tax_table(RA_merged_physeq_differ_f))=="FJ957494.1.1454")
sub_beta_RA_data <- data.table(psmelt(physeq_hiranonis))
sub_beta_RA_data$Timepoint2 = factor(sub_beta_RA_data$Timepoint2, levels =c("Day0", "Day14","Day28","Day42"))
#Figure 6F1
my_comp =list(c("Day0","Day14"),c("Day0","Day28"),c("Day0","Day42"))
ggplot(sub_beta_RA_data, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =0.6) +stat_compare_means(label="p.signif", comparisons=my_comp)+
ylab("Relative Abundance\nof C. hiranonis")+theme_bw()+scale_fill_lancet()+
theme(axis.title.x=element_blank(), legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,0.15), breaks=seq(0,0.1,by=0.02),trans=squish_trans(0,0.002,0.02))


#Hiranonis in 'Complex'
RA_merged_physeq_differ_f=subset_samples(RA_merged_physeq_differ, Remission %in% c("Complex"))
physeq_hiranonis <- subset_taxa(RA_merged_physeq_differ_f, rownames(tax_table(RA_merged_physeq_differ_f))=="FJ957494.1.1454")
sub_beta_RA_data <- data.table(psmelt(physeq_hiranonis))
sub_beta_RA_data$Timepoint2 = factor(sub_beta_RA_data$Timepoint2, levels =c("Day0", "Day14","Day28","Day42"))
#Figure 6F2
my_comp =list(c("Day0","Day14"),c("Day14","Day28"),c("Day28","Day42"))
ggplot(sub_beta_RA_data, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =1) +stat_compare_means(label="p.signif", comparisons=my_comp)+
ylab("Relative Abundance\nof C. hiranonis")+scale_fill_lancet()+theme_bw()+
theme(axis.title.x=element_blank(), legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,0.15), breaks=seq(0,0.1,by=0.02),trans=squish_trans(0,0.002,0.02))



#E. coli in 'Complex'
RA_merged_physeq_differ_f_NDR=subset_samples(RA_merged_physeq_differ, Remission %in% c("Complex","Healthy"))
physeq_Ecoli_NDR <- subset_taxa(RA_merged_physeq_differ_f_NDR, rownames(tax_table(RA_merged_physeq_differ_f_NDR))=="FJ950694.1.1472")
sub_beta_RA_data_NDR <- data.table(psmelt(physeq_Ecoli_NDR))
sub_beta_RA_data_NDR$Timepoint2 = factor(sub_beta_RA_data_NDR$Timepoint2, levels =c("Day0", "Day14","Day28","Day42","Healthy"))
#need paired test
my_comparisons = list (c("Day0","Day14"), c("Day0","Day28"), c("Day14","Day28"), c("Day28","Day42"),c("Day42","Healthy"))
#Figure 4E
ggplot(sub_beta_RA_data_NDR, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =1) +stat_compare_means(comparisons = my_comparisons,label="p.signif")+
ylab("Relative Abundance\nof E. coli")+scale_fill_lancet()+theme_bw()+
theme(axis.title.x=element_blank(), legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,1.8), breaks=seq(0,1,by=0.2),trans=squish_trans(0,0.2,0.1))


#C. perfringens in 'Complex'
physeq_clostri_NDR <- subset_taxa(RA_merged_physeq_differ_f_NDR, rownames(tax_table(RA_merged_physeq_differ_f_NDR))=="New.ReferenceOTU52")
sub_beta_RA_data_NDR <- data.table(psmelt(physeq_clostri_NDR))
sub_beta_RA_data_NDR$Timepoint2 = factor(sub_beta_RA_data_NDR$Timepoint2, levels =c("Day0", "Day14","Day28","Day42","Healthy"))
ggplot(sub_beta_RA_data_NDR, aes(x=Timepoint2, y=Abundance, fill=Timepoint2))+ geom_boxplot(outlier.size =1) +ylab("Relative Abundance\nof Clostridium sp.")+stat_compare_means(comparisons = my_comparisons,label="p.signif")+
scale_fill_lancet()+theme_bw()+theme(axis.title.x=element_blank(), legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 9))+
scale_y_continuous(limits=c(0,0.07))


#comparsion at phylum level
RelativeAb_physeq_ratio_0= transform_sample_counts(merged_physeq_differ, function(x) x / sum(x))
RelativeAb_physeq_ratio =subset_samples(RelativeAb_physeq_ratio_0, Remission%in%c("Diet","Complex"))
RelativeAb_physeq_ratio_phylum_glom = tax_glom(RelativeAb_physeq_ratio, taxrank="Phylum")
RelativeAb_physeq_ratio_phylum_data =data.table(psmelt(RelativeAb_physeq_ratio_phylum_glom))
ratio_plot_data=subset(RelativeAb_physeq_ratio_phylum_data , Phylum%in%c("Fusobacteria", "Firmicutes","Actinobacteria", "Proteobacteria", "Bacteroidetes"))
ratio_plot_data$Remission = gsub("Complex", "NDR", ratio_plot_data$Remission)
ratio_plot_data$Remission = gsub("Diet", "DR", ratio_plot_data$Remission)
Supplementary_Fig_phylum=ggplot(ratio_plot_data, aes(x=Timepoint2,y=Abundance,fill=Remission))+geom_boxplot(outlier.size =1)+scale_fill_lancet()+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size=8), text = element_text(size = 8))+facet_wrap(~Phylum,scales='free_y')+
theme(axis.title.x = element_blank(),strip.background= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")



#preparedata
metabolomic =subset_samples(rared_merged_physeq_f2, Remission %in% c("Diet","Complex","Healthy"))
#Read ND information,  the value that  is below detection (<1 nmol/g) is replaced with '0.5' (half the limit of detection).
ND_infor = read.table (file="ND_information.txt", sep='\t', header =TRUE)
ND_infor =subset (ND_infor, Description =="DogFecal")
ND_infor =subset (ND_infor, Remission %in% c("Diet","Complex","Healthy"))
metabolomic_data = sample_data(metabolomic)
new_metabolomic_data = subset(ND_infor, SampleDescription %in% c(rownames(metabolomic_data)))

#Bile acids # plotting
all_bile=c("AlphamuricholicAcid","BetamuricholicAcid","ChenodeoxycholicAcid_primary","CholicAcid_primary","DeoxycholicAcid","GammamuricholicAcid",
"GlycochenodeoxycholicAcid","GlycocholicAcid","GlycodeoxycholicAcid","LithocholicAcid","OmegamuricholicAcid", "TaurochenodeoxycholicAcid","TaurocholicAcid","TaurodeoxycholicAcid","TaurolithocholicAcid")
#"Diet"
melted_metabolomic_data<- melt (new_metabolomic_data, id.vars=c("DogName","Timepoint2","Remission"),measure.vars=all_bile,variable.name="Bile")
melted_metabolomic_data = subset (melted_metabolomic_data, Remission %in% c("Diet","Healthy"))
melted_metabolomic_data =unique(melted_metabolomic_data)
#For a paired comparison #Missing values for timpoints from "Jack" and "Bebe".
melted_metabolomic_data_test=subset(melted_metabolomic_data, !DogName%in%c("Jack","Bebe") &Timepoint2 %in%c("Day0","Day14","Day42","Healthy"))

#Figure S4
com =list(c("Day0","Day14"),c("Day0","Day42"))
ggplot(melted_metabolomic_data_test, aes(Timepoint2,as.numeric(value)))+geom_boxplot(aes(fill=Timepoint2),outlier.size =1)+geom_signif(textsize=3, vjust=0.6, margin_top = 0.0001, step_increase=0.1, comparisons = com,map_signif_level=TRUE,test.args=list(paired=TRUE))+
ylab("Concentration (nmol/g)")+ theme_bw()+theme(axis.title.x=element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 0.5), text = element_text(size = 10))+theme(legend.position="none")+facet_wrap(~Bile, scales="free_y")+
theme(legend.position="none")+theme(strip.background = element_blank(),legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())


# wilcoxon test #paired test
result<-data.frame("Days0vs14_fc"=integer(),"Days0vs14_pvalue"=integer(),"Days0vs42_fc"=integer(),"Days0vs42_pvalue"=integer())
for (i in all_bile){
day0_primaryBile = subset(melted_metabolomic_data_test, Timepoint2=="Day0" &Bile==i)
day14_primaryBile = subset(melted_metabolomic_data_test, Timepoint2=="Day14" &Bile==i)
day42_primaryBile = subset(melted_metabolomic_data_test, Timepoint2=="Day42" &Bile==i)
foldchange_day0vs14 =mean(as.numeric(day14_primaryBile$value),na.rm=TRUE)/mean(as.numeric(day0_primaryBile$value),na.rm=TRUE)
foldchange_day0vs42 =mean(as.numeric(day42_primaryBile$value),na.rm=TRUE)/mean(as.numeric(day0_primaryBile$value),na.rm=TRUE)
day0VSday14=wilcox.test(as.numeric(day0_primaryBile$value),as.numeric(day14_primaryBile$value), paired=TRUE)
day0VSday42=wilcox.test(as.numeric(day0_primaryBile$value),as.numeric(day42_primaryBile$value),paired=TRUE)
result[i,]<-c(foldchange_day0vs14, day0VSday14$p.value, foldchange_day0vs42, day0VSday42$p.value)
}
#Table S7
write.table(result,"Bile_diet_differ.xls",sep="\t")

#Plot for selected biles
head(melted_metabolomic_data_test)
melted_metabolomic_data_test2 =subset(melted_metabolomic_data_test, Timepoint2 %in%c("Day0","Day14","Day42","Healthy"))
#for paired test
com1=list(c("Day0","Day14"),c("Day0","Day42"))
com2=list(c("Day42","Healthy"))
#DeoxycholicAcid
melted_metabolomic_data3 =subset(melted_metabolomic_data_test2, Bile=="DeoxycholicAcid")
ggplot(melted_metabolomic_data3, aes(x=Timepoint2,y=as.numeric(value),fill=Timepoint2))+geom_boxplot(outlier.size=0.8, alpha=0.9)+
ylab("Concentration (nmol/g stool)")+scale_fill_lancet()+stat_compare_means(comparisons=com1,label="p.signif", paired=TRUE)+stat_compare_means(comparisons=com2,label="p.signif")+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), text = element_text(size =9),legend.position="none")+
theme(strip.background = element_blank(),legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(), plot.title=element_text(vjust=0,hjust=0,size=9))+
ggtitle("Deoxycholic Acid in DR")+scale_y_continuous(limits=c(NA,5000))


#LithocholicAcid
melted_metabolomic_data4 =subset(melted_metabolomic_data_test2, Bile=="LithocholicAcid")
ggplot(melted_metabolomic_data4, aes(Timepoint2,as.numeric(value),fill=Timepoint2))+geom_boxplot(outlier.size=0.8, alpha=0.9)+
ylab("Concentration (nmol/g stool)")+scale_fill_lancet()+stat_compare_means(comparisons=com1,label="p.signif", paired=TRUE)+stat_compare_means(comparisons=com2,label="p.signif")+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), text = element_text(size =9),legend.position="none")+
theme(strip.background = element_blank(),legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(), plot.title=element_text(vjust=0,hjust=0,size=9))+
ggtitle("Lithocholic Acid in DR")+ scale_y_continuous(limits=c(NA,2000))


#**'Complex'
melted_metabolomic_data<- melt(new_metabolomic_data, id.vars=c("DogName","Timepoint2","Remission"),measure.vars=all_bile,variable.name="Bile")
melted_metabolomic_data = subset (melted_metabolomic_data, Remission%in%c("Complex","Healthy"))
melted_metabolomic_data =unique(melted_metabolomic_data)
#
com1=list(c("Day0","Day14"),c("Day0","Day42"))
com2=list(c("Day42","Healthy"))
#Lithocholic acid
complex_lithocholic=subset(melted_metabolomic_data, Bile=="LithocholicAcid"&Timepoint2%in%c("Day0","Day14","Day42","Healthy"))
ggplot(complex_lithocholic, aes(x=Timepoint2, y =as.numeric(value),fill=Timepoint2))+geom_boxplot(outlier.size=0.8,alpha=0.9)+
ylab("Concentration (nmol/g stool)")+scale_fill_lancet()+stat_compare_means(comparisons=com1,label="p.signif", paired=TRUE)+stat_compare_means(comparisons=com2,label="p.signif")+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), text = element_text(size =9),legend.position="none")+
theme(strip.background = element_blank(),legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(), plot.title=element_text(vjust=0,hjust=0,size=9))+
ggtitle("Lithocholic Acid in NDR")+ scale_y_continuous(limits=c(NA,2000))


#deoxycholic acid
complex_deoxycholic=subset(melted_metabolomic_data, Bile=="DeoxycholicAcid"&Timepoint2%in%c("Day0","Day14","Day42","Healthy"))
ggplot(complex_deoxycholic, aes(x=Timepoint2, y =as.numeric(value),fill=Timepoint2))+geom_boxplot(outlier.size=0.8,alpha=0.9)+
ylab("Concentration (nmol/g stool)")+scale_fill_lancet()+stat_compare_means(comparisons=com1,label="p.signif", paired=TRUE)+stat_compare_means(comparisons=com2,label="p.signif")+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), text = element_text(size =9),legend.position="none")+
theme(strip.background = element_blank(),legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(), plot.title=element_text(vjust=0,hjust=0,size=9))+
ggtitle("Deoxycholic Acid in NDR")+scale_y_continuous(limits=c(NA,4800))


#plot heatmap for correlations between bile and OTUs
#Just show the genus that is significantly associated with bile (FDR<0.01)
sig_corr_bile =c("CholicAcid_primary","DeoxycholicAcid","GlycodeoxycholicAcid","LithocholicAcid","TaurocholicAcid")
sig_corr_taxa=c("Bacteroides","Clostridium.sensu.stricto.1","Corynebacterium.1","Escherichia.Shigella","Fusobacterium", "Megamonas",
"Peptoclostridium","Prevotella.9","Staphylococcus","Sutterella")
#
env.taxa.cor_genus$Env =gsub("\\.", "", env.taxa.cor_genus$Env,perl=TRUE)
corr_table_genus= subset(env.taxa.cor_genus, Taxa %in% sig_corr_taxa & Env %in% sig_corr_bile)
corr_table_genus$Env =gsub("Acid", "", corr_table_genus$Env,perl=TRUE)
library(RColorBrewer)
hm.palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(20)
ggplot(corr_table_genus, aes(x=Env, y=Taxa)) + geom_point(aes(size=-log(AdjPvalue), color=Correlation))+
scale_color_gradientn(colours=hm.palette)+
theme_bw()+theme(axis.text.x = element_text(vjust = 0.5), plot.title=element_text(vjust=0,hjust=0,size=9))+ggtitle("Correlation between genus and BA")+
scale_x_discrete(labels=c("Cholic","Deoxycholic", "Glycodeoxycholic","Lithocholic","Taurocholic"))+
theme(legend.margin=margin(t = 0, unit='cm'), legend.key.height=unit(0.6,"line"),legend.key.width=unit(0.6,"line"), axis.text.y=element_text(face="italic"), axis.title.y=element_blank(),axis.title.x=element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 8))



#***************************************************_END_********************************************************************************
