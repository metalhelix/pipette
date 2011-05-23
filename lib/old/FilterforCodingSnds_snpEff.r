
## 1, snps
## 2. indels
source("S:/Bioinformatics/analysis/Gibson/lia/Gibson1/reanalysis/snpEff/snpEff_summary.r");

setwd("S:/Bioinformatics/analysis/Gibson/lia/Gibson1/reanalysis/snpEff")

 
###################  SNPS
outdir<-"S:/Bioinformatics/analysis/Gibson/lia/Gibson1/reanalysis/snpEffResults/"

QL20_snps_ud0 <- read.delim("QL20.snps.snpeff.ud0.collapse.txt")  
QL20_snps_ud0_summary<-snpEff_summary(QL20_snps_ud0, sampleName="QL20_snps",outdir=outdir )

QL6_snps_ud0 <- read.delim("QL6.snps.snpeff.ud0.collapse.txt")  
QL6_snps_ud0_summary<-snpEff_summary(QL6_snps_ud0, sampleName="QL6_snps",outdir=outdir )

QL31_snps_ud0 <- read.delim("QL31.snps.snpeff.ud0.collapse.txt")  
QL31_snps_ud0_summary<-snpEff_summary(QL31_snps_ud0, sampleName="QL31_snps",outdir=outdir )

ICS_wt_snps_ud0 <- read.delim("ICS_wt.snps.filtered.snpeff.collapse.txt")  
ICS_wt_snps_ud0_summary<-snpEff_summary(ICS_wt_snps_ud0, sampleName="ICS_wt_snps",outdir=outdir )

################################# only snps compring to wt
outdir<-"S:/Bioinformatics/analysis/Gibson/lia/Gibson1/reanalysis/snpResultsICS/"

QL20_snps_list<-paste(QL20_snps_ud0[,1], QL20_snps_ud0[,2],QL20_snps_ud0[,3],QL20_snps_ud0[,4], sep="_")
QL6_snps_list<-paste(QL6_snps_ud0[,1], QL6_snps_ud0[,2], QL6_snps_ud0[,3],QL6_snps_ud0[,4],sep="_")
QL31_snps_list<-paste(QL31_snps_ud0[,1], QL31_snps_ud0[,2],QL31_snps_ud0[,3],QL31_snps_ud0[,4] ,sep="_")
wt_snps_list<-paste(ICS_wt_snps_ud0[,1], ICS_wt_snps_ud0[,2], ICS_wt_snps_ud0[,3],ICS_wt_snps_ud0[,4],sep="_")

QL20_snps_ud0<-QL20_snps_ud0[ which( !QL20_snps_list %in% wt_snps_list),]
QL6_snps_ud0<-QL6_snps_ud0[ which( !QL6_snps_list %in% wt_snps_list),]
QL31_snps_ud0<-QL31_snps_ud0[ which( !QL31_snps_list %in% wt_snps_list),]

QL20_snps_ud0_summary<-snpEff_summary(QL20_snps_ud0, sampleName="QL20_snps",outdir=outdir )
QL6_snps_ud0_summary<-snpEff_summary(QL6_snps_ud0, sampleName="QL6_snps",outdir=outdir ) 
QL31_snps_ud0_summary<-snpEff_summary(QL31_snps_ud0, sampleName="QL31_snps",outdir=outdir )

###################  Indels
QL20_indels_ud0 <- read.delim("QL20.indels.snpeff.ud0.collapse.txt")  
QL20_indels_ud0_summary<-snpEff_summary(QL20_indels_ud0, sampleName="QL20_indels",outdir=outdir )

QL6_indels_ud0 <- read.delim("QL6.indels.snpeff.ud0.collapse.txt")  
QL6_indels_ud0_summary<-snpEff_summary(QL6_indels_ud0, sampleName="QL6_indels",outdir=outdir )

QL31_indels_ud0 <- read.delim("QL31.indels.snpeff.ud0.collapse.txt")  
QL31_indels_ud0_summary<-snpEff_summary(QL31_indels_ud0, sampleName="QL31_indels",outdir=outdir )

ICS_wt_indels_ud0 <- read.delim("ICS_wt.indels.filtered.snpeff.collapse.txt")  
ICS_wt_indels_ud0_summary<-snpEff_summary(ICS_wt_indels_ud0, sampleName="ICS_wt_indels",outdir=outdir )


  
QL20_indels_list<-paste(QL20_indels_ud0[,1], QL20_indels_ud0[,2],QL20_indels_ud0[,3],QL20_indels_ud0[,4], sep="_")
QL6_indels_list<-paste(QL6_indels_ud0[,1], QL6_indels_ud0[,2], QL6_indels_ud0[,3],QL6_indels_ud0[,4], sep="_")
QL31_indels_list<-paste(QL31_indels_ud0[,1], QL31_indels_ud0[,2],QL31_indels_ud0[,3],QL31_indels_ud0[,4] ,sep="_")
wt_indels_list<-paste(ICS_wt_indels_ud0[,1], ICS_wt_indels_ud0[,2], ICS_wt_indels_ud0[,3],ICS_wt_indels_ud0[,4],sep="_")

QL20_indels_ud0<-QL20_indels_ud0[ which( !QL20_indels_list %in% wt_indels_list),]
QL6_indels_ud0<-QL6_indels_ud0[ which( !QL6_indels_list %in% wt_indels_list),]
QL31_indels_ud0<-QL31_indels_ud0[ which( !QL31_indels_list %in% wt_indels_list),]

QL20_indels_ud0_summary<-snpEff_summary(QL20_indels_ud0, sampleName="QL20_indels",outdir=outdir )
QL6_indels_ud0_summary<-snpEff_summary(QL6_indels_ud0, sampleName="QL6_indels",outdir=outdir ) 
QL31_indels_ud0_summary<-snpEff_summary(QL31_indels_ud0, sampleName="QL31_indels",outdir=outdir )


