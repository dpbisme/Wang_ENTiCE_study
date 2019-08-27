# Overview of this repo

> This repo houses all of the R code, metadata, and ancillary files for the manuscript "Diet-induced remission in chronic enteropathy is associated with altered microbial community structure and synthesis of secondary bile acids". The 16S rRNA marker gene sequencing data as well as shotgun metagenomic data can be found **[here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA515316/)**. The abstract for the manuscript is copied below:

## Abstract
### Background
The microbiome has been implicated in the initiation and persistence of inflammatory bowel disease. Despite the fact that diet is one of the most potent modulators of microbiome composition and function and that dietary intervention is the first-line therapy for treating pediatric Crohn’s disease, the relationships between diet-induced remission, enteropathy, and microbiome are poorly understood. Here, we leverage a naturally-occurring canine model of chronic inflammatory enteropathy that exhibits robust remission following nutritional therapy, to perform a longitudinal study that integrates clinical monitoring, 16S rRNA gene amplicon sequencing, metagenomic sequencing, metabolomic profiling, and whole genome sequencing to investigate the relationship between therapeutic diet, microbiome, and disease.

### Results
We show that remission induced by a hydrolyzed protein diet is accompanied by alterations in microbial community structure marked by decreased abundance of pathobionts (e.g., *Escherichia coli* and *Clostridium perfringens*), reduced severity of dysbiosis, and increased levels of the secondary bile acids, lithocholic and deoxycholic acid. Physiologic levels of these bile acids inhibited the growth of *E. coli* and *C. perfringens* isolates, *in vitro*. Metagenomic analysis and whole genome sequencing identified the bile acid producer *Clostridium hiranonis* as elevated after dietary therapy and a likely source of secondary bile acids during remission. When *C. hiranonis* was administered to mice, levels of deoxycholic acid were preserved and pathology associated with DSS colitis was ameliorated. Finally, a closely related bile acid producer, *Clostridium scindens*, was associated with diet-induced remission in human pediatric Crohn’s disease. 

### Conclusions 
These data highlight that remission induced by a hydrolyzed protein diet is associated with improved microbiota structure, an expansion of bile acid-producing clostridia, and increased levels of secondary bile acids. Our observations from clinical studies of exclusive enteral nutrition in human Crohn’s disease, along with our in vitro inhibition assays and in vivo studies in mice, suggest that this may be a conserved response to diet therapy with the potential to ameliorate disease. These findings provide insight into diet-induced remission of gastrointestinal disease and could help guide the rational design of more effective therapeutic diets.


> A brief description for each of the files in this repo is included below:

- *ENTiCE_Microbiome_Analysis.R* – All R code used for the analyses reported in the paper
- *ENTiCE_Mapping_file.txt* - mapping file
- *otu_table_mc2_w_tax_no_pynast_failures.json.biom* - OTU table
- *Phylogenetic.fastree* - Phylogenetic tree
- *Bile_acids_MW.txt* - Molecular weight of bile acids
- *Clinical_scores.txt* - CCECAI disease severity scores at start of ENTiCE study (day 0)
- *comparison_of_diet_complex_response.txt* - Genera with differential abundances between diet-responsive (DR) and non-diet-responsive (NDR) dogs.
- *Metabolomics_information.txt* - Fecal metabolite data from ENTiCE study (expressed in nmol/g stool)

```
