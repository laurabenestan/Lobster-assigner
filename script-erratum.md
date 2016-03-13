---
title: "ASSIGNER tutorial for Benestan et al. 's Erratum"
author: "Laura Benestan and Thierry Gosselin"
date: "February 10, 2016"
output:
  html_document:
  theme: united
---
  
  This manual is to help people using **ASSIGNER** to re-analyse the RAD sequencing data from the Benestan et al. paper (2015) published in Molecular Ecology: 
  
  *RAD genotyping reveals fine-scale genetic structuring and provides powerful population assignment in a widely distributed marine species, the American lobster (Homarus americanus).*

This article (Benestan et al. 2015) was published in Molecular Ecology in 2015; since publication, it was brought to our attention that the assignment analysis conducted at the sampling location level (Fig.5) was not reproducible (Eric C. Anderson, personal communication). **Unfortunately, the Training, Holdout, Leave-one-out (THL) method we used for the assignment analysis and described in (Anderson 2010) was not applied correctly at the sampling location level.** 
We mistakenly ranked the markers and conducted the assignment analysis, at the sampling location level, using all the individuals. 

To avoid high-grading bias, the assignment analysis using the THL method requires the samples to be partitioned in a training and holdout datasets (Anderson 2010). Ranking of markers is then conducted on the training samples and the assignment with the selected markers, on the holdout samples. The corrected analysis is shown in this Erratum. The remaining assignment analysis in (Benestan et al. 2015), conducted at the regional levels, were based on the proper THL method.

This Erratum prompted Thierry Gosselin to develop a robust approach in the R language (R Core Team 2015) to achieve reproducible assignment analysis while accounting for the intrinsic nature of GBS/RADseq data analysis (Gosselin et al. 2016). Our approach takes advantage of GSI_SIM, a tool written in C and developed by Eric C. Anderson (Anderson et al. 2008; Anderson 2010) and is implemented in the R package **ASSIGNER** (Gosselin et al. 2016).

The results presented in **this Erratum can be reproducing by following this tutorial**.

Note that if you want to reproduce these results you will want to be on a Mac or a Unix system. You will need to have gsi_sim compiled and on your $PATH.

### *Download the libraries*
```
library(devtools) 
install_github("thierrygosselin/assigner") # to install
install_gsi_sim(fromSource = TRUE) # to install gsi_sim from source
library(assigner)
library(reshape2)
library(ggplot2)
library(stringr)
library(stringi)
library(plyr)
library(dplyr) # load this package after plyr to work properly
library(tidyr)
library(readr)
library(adegenet)
library(randomForestSRC)
library(doParallel)
library(foreach)
library(purrr)
library(utils)
library(iterators)
```

### *Preparing lobster data*

You have to download the file [10156-586.recode.vcf] (http://datadryad.org/resource/doi:10.5061/dryad.q771r) from the Benestan et al. paper's Dryad site cand put it into your current directory within Rstudio.

## *Assessing the extent of high-grading bias*

### *Performing assignment test following THL method for the 11 populations*

We first started by **performing assignment test on the 11 populations** detected by Benestan et al.(2015).

1) Dataset: 10156 SNPs
With the function **GBS_assignment** of the package *ASSIGNER v.0.1.3* we used the filtered the Variant Call Format file (VCF) output by STACKS, which contains only the 10156 markers previously selected in Benestan et al. 2015: `vcf.file = "10156-586.recode.vcf"`. 
Therefore, there is no need to specify a whitelist of markers and a blacklist of individuals then `whitelist.markers = NULL` and `blacklist.id = NULL`. 

2) Keep all the SNPs
Additionally, we can requested if we want to keep only the markers present in all populations. Here, in order to reproduce the same results than the paper **we kept all the 10156 markers**, even if some markers are not missing in some populations. To do so, we specified `common.markers = FALSE`. 

3) Do not correct for LD
Some of the resulting SNPs are positioned on the same 90 pb read. Here, to make the approach reproducible in a pipeline, we evaluated linkage disequilibrium impact on the assignment analysis by using three LD grouping analysis: i) all the markers along the read, independent of their position, ii) the first marker on the read and iii) one marker randomly selected along the read. However, as for Benestan et al paper, **all SNPs were kept even if there were positionned on the same loci**, we specified `snp.ld = NULL`.

4) Use THL method with 50% of the individuals
To reproduce the analysis of Benestan et al. 2015,  we **performed the population assignment test with a THL threshold of 50 %** using the argument `THL = 0.5`. The samples were randomly partitioned, without replacement, into training and holdout samples. The training datasets markers were ranked based on their decreasing global FST. 

5) Subsample population with n = 30
**Since assignment test can be strongly influenced by sample size, we subsample the 11 populations by keeping only 30 individuals per population** (n = 30 being the minimum of individuals present in one population). Then, we specified the argument `subsample = 30`.

6) Do all the subsampling of individuals (required for THL method) 10 times
By subsampling, there is a possibility to create some sampling error (some individuals being samples may be migrants or not). Then **we performed the subsampling 10 times**, which requires the argument:`iterations.subsample = 10`. 

6) Specify the population level used for the assignment test
Here, **we specified the 11 populations** by using the argument `pop.labels = c("TRI", "BON", "GSL", "CAR", "GSL", "GSL", "GSL", "BRA", "GSL", "CAN", "SNS", "SNS", "SNS", "SEA", "BOO", "CCO", "CCO", "RHO")`.
As the name of our population started at the first position and finish at the third position (our population are encoded with three letters, e.g. "TRI") we used: `pop.id.start = 1, pop.id.end = 3`.

7) Do not fill the missing data
No imputations were done on this dataset: `imputations = FALSE`.

8) Run the following command in R 
At the end, we run the command: 

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(data = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.5, blacklist.id = NULL, subsample = 30, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("TRI", "BON", "GSL", "CAR", "GSL", "GSL", "GSL", "BRA", "GSL", "CAN", "SNS", "SNS", "SNS", "SEA", "BOO", "CCO", "CCO", "RHO"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_pop"))
```

### *Performing assignment test with THL 0.5 for the two regions*

We then **performed assignment test on the two regions (North/South)** detected by Benestan et al. (2015). 

1) Specify the region level used for the assignment test 
We run a similar command than the one used for the 11 populations except that here, our population labels corresponds to the two regions (NOR and SOU): `pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU")`. 

2) Do not subsample
Assignment test is not influenced by sample size here since the two regions have a large and similar number of individuals (306 for the North and 280 for the South). Then, we run the argument `subsample = NULL`.

8) Run the following command in R :
Following these arguments, we run the command: 
```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(data = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.5, blacklist.id = NULL, subsample = NULL, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_regions"))
```

### *Performing assignment test with THL 1 for the 11 populations*

We also considered the possibility that dividing the samples into a training and a hold-out or test sets (on average n = 15 each) may have resulted in biased FST values when ranking the markers due to sampling errors associated with too low sample size.

Then **we applied a Leave-one-out (LOO) procedure described in Anderson’s (2010)**, which requires that each individual, in turn, be left out, while the entire process of locus selection and allele frequency estimation is carried out without that individual, and then that individual is assigned back to a population. 

1) Change the method for performing assignment test: use LOO method 
This procedure could be done with the argument `thl = 1`.  

2) Run the following command in R 
For that purpose, we run the following command: 

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(vcf.file = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=NULL, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 1, blacklist.id = NULL, subsample = 30, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("TRI", "BON", "GSL", "CAR", "GSL", "GSL", "GSL", "BRA", "GSL", "CAN", "SNS", "SNS", "SNS", "SEA", "BOO", "CCO", "CCO", "RHO"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_1_populations"))
```

# *Delineating the impact of the number of samples on assignment success*

### *Performing assignment test with 20 individuals subsamples per region

Since regional assignment success was still high after correcting for high-grading bias and samples size were high for both regions, **we used this system (north and south regions) as a positive control to test the impact of the sample size on assignment success**. 

In this case, we expected that regional assignment success, corrected for high-grading bias and performed on the same number of individuals that we sampled per location (which was at maximum 36 individuals), would be similar from the one observed using more individuals. 

1) Use THL method with 50% of the individuals
To test this hypothesis and delineate the gradual effect of the number of samples on assignment success, we selected randomly from 20 to 100 samples from each region, ranked SNPs based on half of these samples (training set), and then calculated the assignment success on the hold-out set (THL method):`THL = 0.50`. 

2) Subsample individuals 10 times
By subsampling, there is a possibility to create some sampling error (some individuals being sampled may be migrants or not). Then **we performed the subsampling 10 times**, which requires the argument:`iterations.subsample = 10`. 

3) Make subsamples of 20,36,50,100 indivioduals
For 20 individuals subsampled, we included the argument `subsample = 20`:

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(data = "10156-586.recode.vcf", whitelist.markers = NULL, snp.LD = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.50, blacklist.id = NULL, subsample = 20, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_regions-sampling20_subsample10"))
```

For 36 individuals subsampled, we included the argument `subsample = 36`:

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(vcf.file = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.5, blacklist.id = NULL, subsample = 36, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_regions-sampling36_subsample10"))
```

For 50 individuals subsampled, we included the argument `subsample = 50`:

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(vcf.file = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.5, blacklist.id = NULL, subsample = 50, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO"), pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_regions-sampling50_subsample10"))
```

For 100 individuals subsampled, we included the argument `subsample = 100`:

```
system.time(assignment.lobster.common.markers.THL.1.sites <- assignment_ngs(vcf.file = "10156-586.recode.vcf", whitelist.markers = NULL, snp.ld = NULL, iterations.subsample=10, common.markers = FALSE, marker.number = c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, "all"), sampling.method = "ranked", thl = 0.5, blacklist.id = NULL, subsample = 50, gsi_sim.filename = "lobster.gsi_sim.txt", keep.gsi.files = FALSE, pop.levels = c("TRI", "BON", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO", "ALM", "ANT", "SJI", "TON"), pop.labels = c("NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "NOR", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "SOU", "ALM", "ANT", "SJI", "TON"), pop.id.start = 1, pop.id.end = 3, imputations = FALSE, parallel.core = 24, folder = "THL_0.5_regions-sampling100_subsample10"))
```

# *Visualisation of the results*

### Figure 1. Proportion of assignment success in relation to the number of markers used and ranked following the THL method (50% of the individuals were used to rank the SNPs) considering 11 populations that were defined in the study. 

We used the output obtained from the population assignment (`THL = 0.50, subsampling = 30`) named `"assignment.results.tsv"`.

First download data and rearrange it:
```
data <- read.table("assignment.results.tsv", header=TRUE, sep="\t",stringsAsFactors=F)
data$MARKER_NUMBER <- as.character(data$MARKER_NUMBER)
neworder <- c("BON","BOO","BRA","CAN","CAR","CCO","GSL","RHO","SEA","SNS","TRI","OVERALL")
data <- arrange(transform(data,
                           CURRENT=factor(CURRENT,levels=neworder)),CURRENT)
```

Then produce the ggplot graph:

```
plot <- ggplot(data, aes(x=MARKER_NUMBER,y=MEAN))+
                 geom_boxplot(aes(fill=CURRENT))+
                 stat_summary(fun.y=mean, geom="point", shape=5, size=3,color="black")+
                 scale_x_discrete(limits=c("100","200","500","1000","2000","3000","4000","5000","6000","7000","8000","9000","10156"))
               x_title="Number of samples"
               y_title="Assignment success (%)"
plot + facet_grid(~CURRENT)+
  facet_wrap(~CURRENT, nrow=4,ncol=3)+
  labs(x=x_title)+
  labs(y=y_title)+
  guides(fill= FALSE, size= FALSE)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="black", linetype="dashed"),
        axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=14, family="Helvetica",face="bold"),
        axis.text.y=element_text(size=14,family="Helvetica",face="bold"))
ggsave("Assignment_pop_THL_col.pdf",width=20,height=25,dpi=300,units="cm",useDingbats=F)
dev.off()
```

### Figure 2. Assignment success in relation to the number of individuals sampled considering the North (panel on the left) and the South (panel in the center) as reference populations.  

We compiled the data obtained from the regional assignment into a file named `"n_effect_summary_regions.txt"`.

First download this file and rearrange data order:
```
data <- read.table("n_effect_summary_regions.txt", header=TRUE)
data$CURRENT <- factor(data$CURRENT,levels = c("NOR", "SOU", "OVERALL"))
neworder <- c("NOR", "SOU", "OVERALL")
data2 <- arrange(transform(data,
                           CURRENT=factor(CURRENT,levels=neworder)),CURRENT)
```

Then use **GGPLOT** package to produce the figure:
```
plot <- ggplot(data, aes(x=NUMBER,y=MEAN, fill=factor(MARKER_NUMBER))+
geom_boxplot()+
stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
scale_x_discrete(limits=c(20,36,50,100,"ALL"))
x_title="Number of samples"
y_title="Assignment success (%)"
```

Then produce a faced grid for each region:
```
plot + facet_grid(~CURRENT)+
  labs(x=x_title)+
  labs(y=y_title)+
  guides(fill= FALSE, size= FALSE)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="black", linetype="dashed"),
        axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=14, family="Helvetica",face="bold"),
        axis.text.y=element_text(size=14,family="Helvetica",face="bold"))
ggsave("Assignment_samples_black_white.pdf",width=20,height=20,dpi=300,units="cm",useDingbats=F)
dev.off()
````
