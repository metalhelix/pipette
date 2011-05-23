

########## take the snpeff.collapse file, classify snps based on it's effect 

snpEff_summary<-function(mydat, sampleName, outdir=getwd()){
   
  ntotal <- nrow(mydat)
  intergenic_snps<-mydat[ (grep ("INTERGENIC", mydat[, "Effect"])),] 
  intronic_snps<-mydat[ (grep ("INTRON", mydat[, "Effect"])),] 

  n_intergenic_snps<-nrow(intergenic_snps)
  n_intronic_snps<-nrow(intronic_snps)  

  exonic_snps<-mydat[ - (grep ("INTERGENIC", mydat[, "Effect"])),] 
  exonic_snps<-exonic_snps[ -(grep ("INTRON", exonic_snps[, "Effect"])),]  ## 
  n_exonic_snps<-nrow(exonic_snps)

  utr_snps <-exonic_snps[ grep("UTR", exonic_snps$Effect),]
  n_utr_snps<-nrow(utr_snps)

  exonic_no_utr_snps<-exonic_snps[ - grep("UTR", exonic_snps$Effect),]

  effect_list<-strsplit(as.character(exonic_no_utr_snps$Effect), split=";")
  ntrans<-unlist(lapply(effect_list, length))
  effect_list_collapse<-data.frame(effect=unlist(effect_list), chrom=rep(exonic_no_utr_snps[,1],ntrans)) 
  
  ids<-unlist(lapply(effect_list, function(x) sum(x %in% c("SYNONYMOUS_CODING", "SYNONYMOUS_START", "SYNONYMOUS_STOP"))))
  nonsyn_snps<-exonic_no_utr_snps[which(ids == 0),]
  
  syn_snps<-exonic_no_utr_snps[which(ids > 0),]
  
  n_nonsyn_snps<-nrow(nonsyn_snps)    
  n_syn_snps<-nrow(syn_snps)   

  more_collapse<-strsplit(as.character(effect_list_collapse[,1]), split=":")
  effect_list_collapse[,1]<-unlist(lapply(more_collapse, function(x) x[1]))
  exonic_effect_summary<-as.data.frame.matrix(table(effect_list_collapse[,1],effect_list_collapse[,2]))
  exonic_effect_summary<-data.frame(type=rownames(exonic_effect_summary), as.data.frame(exonic_effect_summary))

  effect_summary_by_chrom<-table(intergenic_snps$Chrom)
  effect_summary_by_chrom<-rbind(effect_summary_by_chrom, table(intronic_snps$Chrom))
  effect_summary_by_chrom<-rbind(effect_summary_by_chrom, table(utr_snps$Chrom))
  effect_summary_by_chrom<-rbind(effect_summary_by_chrom, table(exonic_no_utr_snps$Chrom))
  effect_summary_by_chrom<-rbind(effect_summary_by_chrom, table(syn_snps$Chrom))
  effect_summary_by_chrom<-rbind(effect_summary_by_chrom, table(nonsyn_snps$Chrom))  
  effect_summary_by_chrom<-data.frame( type=c("intergenic", "intronic", "UTR", "exonic", "syn", "nonsyn"), effect_summary_by_chrom)
  

  write.table(exonic_effect_summary, file=paste(outdir, sampleName, "_exonic_effect_summary.txt", sep=""),quote=F, sep="\t", row.name=F)
  write.table(effect_summary_by_chrom, file=paste(outdir, sampleName, "_effect_summary_by_chrom.txt", sep=""), quote=F, sep="\t", row.name=F)
  write.table(nonsyn_snps, file=paste(outdir, sampleName, "_nonsyn_snps.txt", sep=""), quote=F, sep="\t", row.name=F)
  write.table(utr_snps, file=paste(outdir, sampleName, "_utr_snps.txt", sep=""), quote=F, sep="\t", row.name=F)
  write.table(syn_snps, file=paste(outdir, sampleName, "_syn_snps.txt", sep=""), quote=F, sep="\t", row.name=F)
  write.table(intergenic_snps, file=paste(outdir, sampleName, "_intergenic_snps.txt", sep=""), quote=F, sep="\t", row.name=F)
  write.table(intronic_snps, file=paste(outdir, sampleName, "_intronic_snps.txt", sep=""), quote=F, sep="\t", row.name=F)

  out<-NULL
  out$summary<-effect_summary_by_chrom
  out$exonic_snps_summary<-exonic_effect_summary
  out$nonsyn_snps<-nonsyn_snps
  out$syn_snps<-syn_snps  
  out$utr_snps<-utr_snps
  out$intergenic_snps<-intergenic_snps 
  out$intronic_snps<-intronic_snps
  out
}