#!/usr/bin/Rscript

refgenes = read.table(snakemake@input[['refGene']], sep = "\t", as.is = T, quote = "", header = F)
dim(refgenes)

chrs = refgenes[,3]  # chromosomes, but there is a lot of junk in the chromsome names here that needs to be filtered.

chrs = sapply(strsplit(chrs,"[r_]"),"[", 2)
omits = c(which(chrs == "Y"), which(chrs=="X")) ## drop sex chromosomes

out = cbind(chrs[-omits], refgenes[-omits, 5], refgenes[-omits, 6], refgenes[-omits, 13], refgenes[-omits, 13])

write.table(out, file = snakemake@output[['geneList']], sep = "\t", quote =F, row.name = F, col.name = F)

