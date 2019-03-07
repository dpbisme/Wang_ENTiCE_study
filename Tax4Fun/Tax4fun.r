#!/usr/bin/Rscript
#*****************Set working directory, load packages, import inputfiles, evaluation.
setwd("working_path")
getwd()
#Load libraries
library(biomformat)
library("Tax4Fun")
library("phyloseq")
library("ggsci")
library(dplyr)
library("tibble")
library("ggpubr")

#Import mapping file
mapping <- import_qiime_sample_data(mapfilename = "MappingFile_tax4Fun.txt")
qiimebio<-importQIIMEBiomData("Tax4Fun_otu_table_collapsed.json.biom")
referencedata <- "E:\\นคื๗\\Upenn-2017-4-18\\experiments\\canine_IBD\\MS\\2019-2-28-code_sub\\Tax4Fun\\SILVA123\\SILVA123"

#Kegg orthologs (KO) analysis
Tax4FunOutput <- Tax4Fun(qiimebio, referencedata)
print(Tax4FunOutput)
Tax4FunProfile <- data.frame(t(Tax4FunOutput$Tax4FunProfile))
#save table
write.table(Tax4FunProfile,"Tax4FunProfile_Export.xls",sep="\t")

# KEGG Pathway analysis
KEGG_pathways <- Tax4Fun(qiimebio, fctProfiling = FALSE, refProfile = "PAUDA", normCopyNo = TRUE, folderReferenceData = referencedata)
print(KEGG_pathways)
#
rownames(KEGG_pathways$Tax4FunProfile)
KEGG_pathwaysProfile <- data.frame(t(KEGG_pathways$Tax4FunProfile))
#save table
write.table(KEGG_pathwaysProfile,"KEGG_pathwaysProfile_Export.xls",sep="\t")

#Get the mapping data only for diet responsive dog
diet_samples0=subset(mapping, Remission %in% c("Diet","Healthy"))
diet_samples = subset (diet_samples0, Timepoint2 %in% c("Day0", "Day14","Day28","Day42","Healthy"))
#get Tax4funProfile data only for diet-responsive dog
Tax4FunProfile_diet=Tax4FunProfile[,rownames(diet_samples)]
#get only data for fecal samples
Fecal_sample_diet <- grep("DogFecal",colnames(Tax4FunProfile_diet),perl=TRUE,value=TRUE)
#get individual data for days
day0 <- grep("_d0_",Fecal_sample_diet,perl=TRUE, value=TRUE)
day14 <- grep("_d14_",Fecal_sample_diet,perl=TRUE,value=TRUE)
day28 <- grep("_d28_",Fecal_sample_diet,perl=TRUE,value=TRUE)
day42<- grep("_d42_",Fecal_sample_diet,perl=TRUE,value=TRUE)

#For a paired test, you need to check the orders and elements on the datasets you want to compare
day0
day14
#Set up data frame of result
result<-data.frame("Days0vs14_lfc"=integer(),"Days0vs14_pvalue"=integer())
#loop for wilcoxon test
for(i in 1:nrow(Tax4FunProfile_diet)){
data_a = as.numeric(c(Tax4FunProfile_diet[i, day0]))
data_a = as.numeric(append (data_a, 'NA', 7))
data_b = as.numeric(c(Tax4FunProfile_diet[i, day14]))
b_a_foldchange_log2 =log2(mean(data_b, na.rm=TRUE)/mean(data_a, na.rm=TRUE))
b_vs_a=wilcox.test(data_b, data_a, alternative = "two.sided",paired=TRUE, exact=FALSE)
result[i,]<-c(b_a_foldchange_log2, b_vs_a$p.value)
}
#
rownames(result)=rownames((Tax4FunProfile_diet))
#Adjusted p-values #BH
result$Days0vs14_fdr =p.adjust(result$Days0vs14_pvalue, method = "BH")
write.table(result[,c(1,2,3)],"Tax4Fun_KO_diff_days0vs14.xls",sep="\t")

#PCA analysis
Tax4FunProfile_diet_PCA <- data.frame(t(Tax4FunProfile_diet))
Tax4FunProfile_diet_PCA<- rownames_to_column(Tax4FunProfile_diet_PCA, 'description')
Tax4FunProfile_diet_fecal_PCA<- subset(Tax4FunProfile_diet_PCA, description %in% Fecal_sample_diet)
Tax4FunProfile_diet_fecal_PCA <- mutate(Tax4FunProfile_diet_fecal_PCA, group= ifelse(description%in% day0,"Day0", ifelse(description %in% day14, "Day14",ifelse(description %in%day28, "Day28", ifelse(description %in%day42, "Day42", "Healthy")))))
Tax4FunProfile_diet_fecal_PCA<-subset(Tax4FunProfile_diet_fecal_PCA, group %in% c("Day0","Day14"))
KO_groups<- Tax4FunProfile_diet_fecal_PCA$group
Tax4FunProfile_diet_fecal_PCA =subset (Tax4FunProfile_diet_fecal_PCA,select =-c(group,description))
Tax4FunProfile_diet_fecal_PCA_f<-Tax4FunProfile_diet_fecal_PCA[, colSums(Tax4FunProfile_diet_fecal_PCA)>0]
#For KO
require(factoextra)
#Not scale for
pca.res <- prcomp(na.omit(Tax4FunProfile_diet_fecal_PCA_f), scale.=F, retx=T)
fviz_pca_ind(pca.res, label="none",axes = c(1, 2), alpha.ind =1,  habillage=KO_groups, invisible="quali",
             addEllipses=FALSE, ellipse.level=0.66)+scale_color_lancet()+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+theme(text = element_text(size = 12))+ggtitle("")

#get Tax4fun data only for diet-responsive dog
KEGG_pathwaysProfile_diet=KEGG_pathwaysProfile[,rownames(diet_samples)]
#get only data for fecal samples
Fecal_sample_diet_KEGG <- grep("DogFecal",colnames(KEGG_pathwaysProfile_diet),perl=TRUE, value =TRUE)
day0 <- grep("_d0_",Fecal_sample_diet_KEGG,perl=TRUE, value=TRUE)
day14 <- grep("_d14_",Fecal_sample_diet_KEGG,perl=TRUE,value=TRUE)
day28 <- grep("_d28_",Fecal_sample_diet_KEGG,perl=TRUE,value=TRUE)
day42<- grep("_d42_",Fecal_sample_diet_KEGG,perl=TRUE,value=TRUE)

#For a paired test, you need to check the orders and elements on the datasets you want to compare
day0
day14
result_pathway<-data.frame("Days0vs14_lfc"=integer(),"Days0vs14_pvalue"=integer())
#loop for wilcoxon test
for(i in 1:nrow(KEGG_pathwaysProfile_diet)){
data_a = as.numeric(c(KEGG_pathwaysProfile_diet[i, day0]))
#for paired comparison
data_a = as.numeric(append (data_a, 'NA', 7))
data_b = as.numeric(c(KEGG_pathwaysProfile_diet[i, day14]))
b_a_foldchange_log2 =log2(mean(data_b,na.rm=TRUE)/mean(data_a,na.rm=TRUE))
b_vs_a=wilcox.test(data_b,data_a, alternative = "two.sided",paired=TRUE)
result_pathway[i,]<-c(b_a_foldchange_log2, b_vs_a$p.value)
}
#
rownames(result_pathway)=rownames((KEGG_pathwaysProfile_diet))
#FDR values added
result_pathway$Days0vs14_fdr =p.adjust(result_pathway$Days0vs14_pvalue, method = "fdr")
write.table(result_pathway[,c(1,2,3)],"Tax4Fun_pathway_diff_days0vs14.xls",sep="\t")

#PCA analysis of KEGG pathways
KEGG_pathwaysProfile_diet_PCA <- data.frame(t(KEGG_pathwaysProfile_diet))
KEGG_pathwaysProfile_diet_PCA<- rownames_to_column(KEGG_pathwaysProfile_diet_PCA, 'description')
KEGG_pathwaysProfile_diet_fecal_PCA<- subset(KEGG_pathwaysProfile_diet_PCA, description %in% Fecal_sample_diet)
KEGG_pathwaysProfile_diet_fecal_PCA <- mutate(KEGG_pathwaysProfile_diet_fecal_PCA, Timepoint= ifelse(description%in% day0,"Day0", ifelse(description %in% day14, "Day14",ifelse(description %in%day28, "Day28", ifelse(description %in%day42, "Day42", "Healthy")))))
KEGG_pathwaysProfile_diet_fecal_PCADays0vs14 <- subset(KEGG_pathwaysProfile_diet_fecal_PCA, Timepoint %in% c("Day0", "Day14"))
KEGG_groups <- KEGG_pathwaysProfile_diet_fecal_PCADays0vs14$Timepoint
KEGG_pathwaysProfile_diet_fecal_PCADays0vs14 =subset (KEGG_pathwaysProfile_diet_fecal_PCADays0vs14,select =-c(Timepoint,description))
KEGG_pathwaysProfile_diet_fecal_PCA_f<-KEGG_pathwaysProfile_diet_fecal_PCADays0vs14[, colSums(KEGG_pathwaysProfile_diet_fecal_PCADays0vs14)>0]
#For KEGG pathway
#Figure 5A
require(factoextra)
#Not scale
pca.res <- prcomp(KEGG_pathwaysProfile_diet_fecal_PCA_f, scale.=F, retx=T)
fviz_pca_ind(pca.res, label="none",axes = c(1, 2), alpha.ind =1,  habillage=KEGG_groups,invisible="quali",pointsize = 1.5, title="PCA of pathways", addEllipses = TRUE,
            ellipse.level=0.66,palette = c("#00AFBB", "#FC4E07"))+theme_bw()+
	     theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
	     theme(text = element_text(size = 9), legend.position=c(0.84,0.85),legend.margin=margin(t = 0, unit='cm'), legend.title=element_blank(),legend.key.width=unit(0.8,"line"),legend.key.height=unit(0.6,"line"))

#Contributions of each dim
fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 100))
fviz_contrib(pca.res, choice = "var", axes = 1:2, top=15)

#plot for dim
KEGG_pathwaysProfile_diet_fecal_PCA =subset(KEGG_pathwaysProfile_diet_fecal_PCA, Timepoint != "Healthy")
Timepoints= KEGG_pathwaysProfile_diet_fecal_PCA$Timepoint
KEGG_pathwaysProfile_diet_fecal_PCA_Days =subset (KEGG_pathwaysProfile_diet_fecal_PCA, select =-c(Timepoint,description))
KEGG_pathwaysProfile_diet_fecal_PCA_Days_f<-KEGG_pathwaysProfile_diet_fecal_PCA_Days[, colSums(KEGG_pathwaysProfile_diet_fecal_PCA_Days)>0]
pca.res_all <- prcomp(KEGG_pathwaysProfile_diet_fecal_PCA_Days_f, scale.=F, retx=T)
pca_all = data.frame(pca.res_all$x, Timepoints)
fviz_contrib(pca.res_all, choice = "var", axes = 1, top=15)
com1=list(c("Day0","Day14"),c("Day0","Day28"), c("Day0","Day42"))
#Figure 5B
ggplot(pca_all, aes(x=Timepoints, y = PC1, fill=Timepoints))+geom_boxplot(outlier.size =0.8)+stat_compare_means(comparisons=com1,label="p.signif")+
theme_bw()+ylab("Dim 1")+
theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), legend.position="none", strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 9),plot.title=element_text(hjust=0,size=9))+
scale_fill_lancet()+guides(fill=FALSE)+ggtitle("Dim 1")+scale_y_continuous(limits=c(NA, 0.08))

#Heatmap of selected pathways
selected_pathways = c("ko00071; Fatty acid metabolism", "ko00100; Steroid biosynthesis", "ko00281; Geraniol degradation", "ko01040; Biosynthesis of unsaturated fatty acids","ko00592; alpha-Linolenic acid metabolism",
"ko04973; Carbohydrate digestion and absorption", "ko00052; Galactose metabolism", "ko00511; Other glycan degradation", "ko00010; Glycolysis / Gluconeogenesis","ko00051; Fructose and mannose metabolism","ko00121; Secondary bile acid biosynthesis")
selected_pathways_diet = KEGG_pathwaysProfile_diet[selected_pathways,]
selected_pathways_diet_fecal = selected_pathways_diet[,grep ("DogFecal", colnames(selected_pathways_diet))]
z.mat <- t(scale(t(selected_pathways_diet_fecal), scale=TRUE, center=TRUE))
selected_pathways_diet_fecal_heatmap = data.frame(t(z.mat), description = colnames(z.mat))
#for order
selected_pathways_diet_fecal_heatmap2 <- plyr::mutate (selected_pathways_diet_fecal_heatmap, Timepoint= ifelse(description%in% day0,"Day0", ifelse(description %in% day14, "Day14",ifelse(description %in%day28, "Day28", ifelse(description %in%day42, "Day42", "Healthy")))))
selected_pathways_diet_fecal_heatmap3 = subset(selected_pathways_diet_fecal_heatmap2, Timepoint %in%c("Day0", "Day14"))
selected_pathways_diet_fecal_heatmap4= data.table::melt(selected_pathways_diet_fecal_heatmap3, id.vars = c("Timepoint","description"),variable.name="KEGG_pathways")
selected_pathways_diet_fecal_heatmap5= selected_pathways_diet_fecal_heatmap4[order(selected_pathways_diet_fecal_heatmap4$Timepoint),]
selected_pathways_diet_fecal_heatmap5$description = factor (selected_pathways_diet_fecal_heatmap5$description, levels = unique(selected_pathways_diet_fecal_heatmap5$description))
selected_pathways_diet_fecal_heatmap5$KEGG_pathways = gsub("(\\S+\\.\\.)", "", selected_pathways_diet_fecal_heatmap5$KEGG_pathways, perl=TRUE)
selected_pathways_diet_fecal_heatmap5$KEGG_pathways = gsub("(\\.)", " ", selected_pathways_diet_fecal_heatmap5$KEGG_pathways, perl=TRUE)
order_pathways = selected_pathways = c("Fatty acid metabolism", "Steroid biosynthesis", "Geraniol degradation", "Biosynthesis of unsaturated fatty acids","alpha Linolenic acid metabolism",
"Carbohydrate digestion and absorption", "Galactose metabolism", "Other glycan degradation", "Gluconeogenesis","Fructose and mannose metabolism","Secondary bile acid biosynthesis")
selected_pathways_diet_fecal_heatmap5$KEGG_pathways =factor(selected_pathways_diet_fecal_heatmap5$KEGG_pathways, levels=order_pathways)
#Figure 5C
library(RColorBrewer)
hm.palette <- colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20)
ggplot(selected_pathways_diet_fecal_heatmap5, aes(x=description, y=KEGG_pathways, fill=value)) + geom_tile() +
      scale_fill_gradientn(colours = hm.palette,limits=c(-3.07, 3),na.value= "#A50026")+
        theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(),axis.text.y=element_text(size=9))+
	theme(text = element_text(size = 9),legend.margin=margin(t = 0, unit='cm'), legend.key.height=unit(1,"line"),legend.key.width=unit(0.6,"line"),plot.title=element_text(vjust=0, hjust=-1, size=9))+ggtitle("KEGG Pathways")
 # coord_equal()

##Secondary bile acid biosynthesis pathway
secondary_BA = "ko00121; Secondary bile acid biosynthesis"
selected_pathways_BA = KEGG_pathwaysProfile_diet[secondary_BA,]
selected_pathways_diet_BA= selected_pathways_BA[,grep ("DogFecal", colnames(selected_pathways_BA))]
selected_pathways_diet_fecal_BA = data.frame(t(selected_pathways_diet_BA), description = colnames(selected_pathways_diet_BA))
#Sort
selected_pathways_diet_fecal_BA2 <- plyr::mutate (selected_pathways_diet_fecal_BA, Timepoint= ifelse(description%in% day0,"Day0", ifelse(description %in% day14, "Day14",ifelse(description %in%day28, "Day28", ifelse(description %in%day42, "Day42", "Healthy")))))
selected_pathways_diet_fecal_BA2 = selected_pathways_diet_fecal_BA2 [order(selected_pathways_diet_fecal_BA2$Timepoint),]
selected_pathways_diet_fecal_BA3  = subset(selected_pathways_diet_fecal_BA2 , Timepoint %in%c("Day0", "Day14","Day28","Day42"))
#Figure 5D
library(ggpubr)
my_comparisons= list (c("Day0", "Day14"), c("Day0", "Day28"), c("Day0","Day42"))
ggplot(selected_pathways_diet_fecal_BA3, aes(x= Timepoint, y = ko00121..Secondary.bile.acid.biosynthesis, fill=Timepoint))+geom_boxplot(outlier.size =0.8)+stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "bottom")+
ylab("Relative Abundance")+theme_bw()+
theme(legend.position="none", strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 9),plot.title=element_text(vjust=0, size=9))+ggtitle("Secondary BA synthesis")+
scale_y_continuous(limits=c(NA,8e-4))+scale_fill_lancet()