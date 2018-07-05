#!/usr/bin/Rscript

### R script for generating setIDs which can be fed to SKAT

sets = scan(snakemake@input[['set']], what = "char")

setname = sets[1]
i = 2

sink(snakemake@output[['setid']])

while (1) {
  while (sets[i] != "END") {
      cat(setname,"\t",sets[i],"\n", sep = "")
      i = i+1
  }
  if (i == length(sets))
      break;
  i = i+1
  setname = sets[i]
  i = i+1
}
sink()
