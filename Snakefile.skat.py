from snakemake.utils import R

### Set constants
BEDFOLDER 	= config['bedFolder']
REFGENE		= config['refGene']
SSDFOLDER 	= 'ssd/'
PSSDFOLDER	= 'pheno_ssd/'
PLINK		= 'plink '
PHENOFILE	= 'gwas.pheno.covar.quantNorm.numeric.txt'#'covars.quantnorm.txt' #'quant.rank.txt'#'covars.log2.txt' #'covarnotransform'#'pheno.final.csv'#validPhenos.txt'
PHENOSTOREMOVE = 'samplesForRaw.txt'#"samplesForWhichPhenosAreGiven.txt"#'phenosgivencovarnotransform' #

PHENOTYPES 	= ['solAB42','insolAB42', 'APPBAPPA', 'sAPPalpha', 'sAPPbeta'] #'AB42AB40', 'AB40', 'AB42']

SKAT_FILES	= ['SSD', 'Info']
BED_FILES 	= ['bed', 'bim', 'fam']

SKATMAF = 0.1

### Set startup configuration
exec(open(config['runConfiguration']).read())


onstart:
	print("Starting pipeline, datasets: " + ','.join(DATASETS) + "\nChromosomes: " + ', '.join([str(c) for c in CHR]))	

onsuccess:
	print("Finished.")

	
rule all:
	input	: expand('assoc/{dataset}.skat', dataset=DATASETS),
			  expand('assoc/{dataset}.crskat', dataset=DATASETS),#, 'assoc/NBB.omni.skat'#, 'assoc/NBB.omni.skat'#expand(SSDFOLDER+ '{dataset}.chr{chr}.{ext}', dataset=DATASETS, chr=CHR, ext=SKAT_FILES) 
			  expand('assoc/{dataset}.{pheno}.skat', dataset=DATASETS, pheno=PHENOTYPES),
			  expand('assoc/{dataset}.{pheno}.crskat', dataset=DATASETS, pheno=PHENOTYPES),
			  'xls/skat.xls', 'xls/crskat.xls',
			  expand('genoPlots/{dataset}.{pheno}.{skat}.pdf', dataset=DATASETS, pheno=PHENOTYPES, skat=['skat','crskat'])
			  #'pheno_chr_skat_extended/NBB.exome.chr9.insolAB42.crskat'
rule clean:
	shell	: 'rm -rf tmp/ ssd/ meta/ assoc/ chr_skat/ pheno_assoc/ pheno_chr_skat/ pheno_ssd/ pheno_bed/'

rule makeXLS:
	input	: expand('assoc/{dataset}.skat', dataset=DATASETS), expand('assoc/{dataset}.{pheno}.skat', dataset=DATASETS, pheno=PHENOTYPES)
	output	: xls='xls/skat.xls'
	script	: 'scripts/makeXLS.R'
	
rule makeXLSCR:
	input	: expand('assoc/{dataset}.crskat', dataset=DATASETS), expand('assoc/{dataset}.{pheno}.crskat', dataset=DATASETS, pheno=PHENOTYPES)
	output	: xls='xls/crskat.xls'
	script	: 'scripts/makeXLS.R'	

rule mergePlots:
	input: 'genoPlots/{dataset}.{pheno}.{skattype}/'
	output: 'genoPlots/{dataset}.{pheno}.{skattype}.pdf'
	shell: 'convert {input}* {output}'	
	
rule plotBest:
	input: assoc='assoc/{dataset}.{pheno}.{skattype}', pheno=PHENOFILE
	output: 'genoPlots/{dataset}.{pheno}.{skattype}/'
	run:
		R("""
		
		library(SKAT)
		library(stringr)
		library(plotrix)
		
		assoc <- read.delim('{input.assoc}', sep=" ")
		
		n <- min(nrow(assoc),16)
		
		if (n > 0) {{
		for (i in 1:n) {{
		
			row <- assoc[i, ]
			ssd = paste('ssd/{wildcards.dataset}.chr',row$CHR,'.SSD', sep="")
			fam = paste('bed/{wildcards.dataset}.chr',row$CHR,'.fam', sep="")
			info = paste('ssd/{wildcards.dataset}.chr',row$CHR,'.Info', sep="")
			pheno <- '{input.pheno}'
			
			FAM<- Read_Plink_FAM_Cov(fam, pheno, Is.binary=FALSE)
			y <- FAM${wildcards.pheno}			
			
			SSD.INFO = Open_SSD(File.SSD=ssd, File.Info=info)

			SetIndex <- which(SSD.INFO$SetInfo$SetID == row$SetID)
			
			genos <- Get_Genotypes_SSD(SSD.INFO, SetIndex, TRUE)
			
			genos <- as.data.frame(genos)
			
			genos${wildcards.pheno} <- y
			
			genos$haplotypes <- do.call(paste0, genos[colnames(genos)[1:length(colnames(genos))-1]])
			
			genos$haplotypes <- as.factor(genos$haplotypes)

			png(paste('{output}',str_pad(i, 2, pad="0"),'_',row$SetID,'.png',sep=""))
			boxplot(genos${wildcards.pheno} ~ genos$haplotypes, main=paste('{wildcards.dataset}',' chr', row$CHR, sep=""), ylab='{wildcards.pheno}', xlab=row$SetID) 
			colors <- ifelse(FAM$braak > 3, 'red', 'green')
			points(genos${wildcards.pheno} ~ genos$haplotypes, col=colors)
			mtext(paste(colnames(genos)[1:max(1, ncol(genos)-2)], collapse=" "),side=3,outer=F) 
			
			
			# count genotypes
			f <- summary(genos$haplotypes)
			f <- as.data.frame(f)
			f <- cbind(rownames(f),f$f)
			colnames(f) <- c('Hap','n')
			f <- apply(f, 2, str_pad, width=nchar(f[1,1]), pad=" ")
			addtable2plot(x="bottomright",table=f)
			
			dev.off()
		}}
		
		}}
		
		
		""")
	
	
rule multipleTesting:
	input: 'assoc_pre/{dataset}.{extension}'
	output: filtered='assoc/{dataset}.{extension}', full='assoc_full/{dataset}.{extension}'
	run:
		R("""
		
		assoc <- read.delim('{input}', sep=" ")
		
		assoc$adj.p <- p.adjust(assoc$P.value, method="fdr")
		assoc <- assoc[order(assoc$P.value),]
		
		write.table(assoc, file='{output.full}', row.names=F, quote=F)
		
		assoc <- assoc[assoc$adj.p <= 0.5,]
		
		assoc <- assoc[order(assoc$adj.p),]
		
		write.table(assoc, file='{output.filtered}', row.names=F, quote=F)
		
		""")
	
rule mergeSKAT:
	input	: expand('chr_skat_extended/{{dataset}}.chr{chr}.skat',chr=CHR)
	output	: 'assoc_pre/{dataset}.skat'
	shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 1) print;}}' | sort -k9 -g > {output}"

rule mergeSKATCR:
	input	: expand('chr_skat_extended/{{dataset}}.chr{chr}.crskat',chr=CHR)
	output	: 'assoc_pre/{dataset}.crskat'
	#shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 0.05) print;}}' | sort -k8 -g > {output}"
	shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 1) print;}}' | sort -k8 -g > {output}"

	
rule mergePhenoSKAT:
	input	: expand('pheno_chr_skat_extended/{{dataset}}.chr{chr}.{{pheno}}.skat', chr=CHR)
	output	: 'assoc_pre/{dataset}.{pheno}.skat'
	#shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 0.05) print;}}' | sort -k5 -g > {output}"
	shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 1) print;}}' | sort -k5 -g > {output}"
	
	
rule mergePhenoSKATCR:
	input	: expand('pheno_chr_skat_extended/{{dataset}}.chr{chr}.{{pheno}}.crskat', chr=CHR)
	output	: 'assoc_pre/{dataset}.{pheno}.crskat'
	#shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 0.05) print;}}' | sort -k8 -g > {output}"
	shell 	: "awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}	' {input} | awk '{{ if (NR==1 || $2 <= 1) print;}}' | sort -k8 -g > {output}"
	
rule appendChr_SKAT:
	input: 'chr_skat/{dataset}.chr{chr}.skat'
	output: 'chr_skat_extended/{dataset}.chr{chr}.skat'
	#shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 0.05) print $0 " {wildcards.chr}" }}' > {output}"""
	shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 1) print $0 " {wildcards.chr}" }}' > {output}"""


rule appendChr_SKATCR:
	input: 'chr_skat/{dataset}.chr{chr}.crskat'
	output: 'chr_skat_extended/{dataset}.chr{chr}.crskat'
	#shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 0.05) print $0 " {wildcards.chr}" }}' > {output}"""
	shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 1) print $0 " {wildcards.chr}" }}' > {output}"""
	

rule appendChr_Pheno_SKAT:
	input: 'pheno_chr_skat/{dataset}.chr{chr}.{pheno}.skat'
	output: 'pheno_chr_skat_extended/{dataset}.chr{chr}.{pheno}.skat'
	#shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 0.05) print $0 " {wildcards.chr}" }}' > {output}"""
	shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 1) print $0 " {wildcards.chr}" }}' > {output}"""

	
rule appendChr_Pheno_SKATCR:
	input: 'pheno_chr_skat/{dataset}.chr{chr}.{pheno}.crskat'
	output: 'pheno_chr_skat_extended/{dataset}.chr{chr}.{pheno}.crskat'
	#shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 0.05) print $0 " {wildcards.chr}" }}' > {output}"""
	shell: """awk 'FNR==1 && NR!=1 {{ while (/^SetID/) getline; }} 1 {{print}}' {input} | awk '{{ if (NR==1) print $0 " CHR"; if ($2 <= 1) print $0 " {wildcards.chr}" }}' > {output}"""
	

rule runSKAT:
	input	: fam	= BEDFOLDER + '{dataset}.chr{chr}.fam', 
			  ssd 	= SSDFOLDER + '{dataset}.chr{chr}.SSD',
			  info	= SSDFOLDER + '{dataset}.chr{chr}.Info',
			  pheno = PHENOFILE
			  
	output	: assoc='chr_skat/{dataset}.chr{chr}.skat',cr='chr_skat/{dataset}.chr{chr}.crskat', skatqqplot='qqplots/skat.{dataset}.chr{chr}.png'#, crskatqqplot='qqplots/crskat.{dataset}.chr{chr}.png'
	params  : skatmaf=SKATMAF
	run		: 
		R("""
		
		library("SKAT")
		
		fam 	<- '{input.fam}'
		ssd  	<- '{input.ssd}'
		info	<- '{input.info}'		
		pheno 	<- '{input.pheno}'
		
		SSD.INFO<- Open_SSD(File.SSD=ssd, File.Info=info)
		FAM		<- Read_Plink_FAM_Cov(fam, pheno)
		

		
		#FILTER SSD INFOS
		
		newSSD.Info = SSD.INFO
		newSSD.Info$SetInfo <- newSSD.Info$SetInfo[which(newSSD.Info$SetInfo$SetSize >= 4),]
		newSSD.Info$nSets <- nrow(newSSD.Info$SetInfo)
		newSSD.Info$nSNPs <- sum(newSSD.Info$SetInfo$SetSize)

		#SSD.INFO = newSSD.Info
		
		
		y 		<- FAM$Phenotype
		
		#age <- FAM$age
		#gender <- FAM$gender
		#braak <- FAM$braak 
		
		obj 	<- SKAT_Null_Model(y ~ age + gender + braak, data=FAM, out_type="D")				
		#assoc 	<- SKAT.SSD.All(SSD.INFO, obj)
		assoc 	<- SKATBinary.SSD.All(SSD.INFO, obj, kernel="linear.weighted", 	method="optimal.adj", max_maf={params.skatmaf})
		common.rare <- SKAT_CommonRare.SSD.All(SSD.INFO, obj, method="optimal.adj")#, max_maf={params.skatmaf})

		assoc$results <- assoc$results[which(assoc$results$N.Marker.Test >= 4),]
		common.rare$results <- common.rare$results[which(common.rare$results$N.Marker.Test >= 4),]
				
		#assoc$results$adj.p <- p.adjust(assoc$results$P.value, method="fdr")
		#common.rare$results$adj.p <- p.adjust(common.rare$results$P.value, method="fdr")

		Get_EffectiveNumberTest(assoc$results$MAP, alpha=0.05)
		#Get_EffectiveNumberTest(common.rare$results$MAP, alpha=0.05)
		
		png(filename="{output.skatqqplot}")
		QQPlot_Adj(assoc$results$P.value, assoc$results$MAP)	
		dev.off()
		
		#png(filename="")
		#QQPlot_Adj(common.rare$results$P.value, common.rare$results$MAP)
		#dev.off()
		
		write.table(assoc$results, file='{output.assoc}', row.names=F, quote=F)
		write.table(common.rare$results, file='{output.cr}', row.names=F, quote=F)
		""")

rule runPhenoSKAT:
	input	: fam	= 'pheno_bed/{dataset}.chr{chr}.fam', 
			  ssd 	= PSSDFOLDER + '{dataset}.chr{chr}.SSD',
			  info	= PSSDFOLDER + '{dataset}.chr{chr}.Info',
			  pheno = PHENOFILE
	output	: assoc='pheno_chr_skat/{dataset}.chr{chr}.{pheno}.skat',cr='pheno_chr_skat/{dataset}.chr{chr}.{pheno}.crskat'#, skatqqplot='qqplots/skat.{dataset}.chr{chr}.{pheno}.png', crskatqqplot='qqplots/crskat.{dataset}.chr{chr}.{pheno}.png'
	params : skatmaf=SKATMAF
	run		: 
		R("""
		
		library("SKAT")
		
		fam 	<- '{input.fam}'
		ssd  	<- '{input.ssd}'
		info	<- '{input.info}'
		pheno 	<- '{input.pheno}'
		
		SSD.INFO<- Open_SSD(File.SSD=ssd, File.Info=info)
		FAM		<- Read_Plink_FAM_Cov(fam, pheno, Is.binary=FALSE)
		
		#FILTER SSD INFOS
		
		newSSD.Info = SSD.INFO
		newSSD.Info$SetInfo <- newSSD.Info$SetInfo[which(newSSD.Info$SetInfo$SetSize >= 4),]
		newSSD.Info$nSets <- nrow(newSSD.Info$SetInfo)
		newSSD.Info$nSNPs <- sum(newSSD.Info$SetInfo$SetSize)

		#SSD.INFO = newSSD.Info
		
		y 		<- FAM${wildcards.pheno}
		#print(head(FAM$insolAB42))
		#print(y)
		
		#age <- FAM$age
		#gender <- FAM$gender
		#braak <- FAM$braak 
		
		obj 	<- SKAT_Null_Model({wildcards.pheno} ~ age + gender + braak, data=FAM, out_type="C")		
		assoc 	<- SKAT.SSD.All(SSD.INFO, obj, method="optimal.adj", max_maf={params.skatmaf})		
		common.rare <- SKAT_CommonRare.SSD.All(SSD.INFO, obj, method="optimal.adj")#, max_maf={params.skatmaf})
		head(common.rare)
		
		assoc$results <- assoc$results[which(assoc$results$N.Marker.Test >= 4),]
		common.rare$results <- common.rare$results[which(common.rare$results$N.Marker.Test >= 4),]
		
		#assoc$results$adj.p <- p.adjust(assoc$results$P.value, method="fdr")
		#common.rare$results$adj.p <- p.adjust(common.rare$results$P.value, method="fdr")

		#Get_EffectiveNumberTest(assoc$results$MAP, alpha=0.05)
		#Get_EffectiveNumberTest(common.rare$results$MAP, alpha=0.05)
		
		#png(filename="")
		#QQPlot_Adj(assoc$results$P.value, assoc$results$MAP)	
		#dev.off()
		
		#png(filename="")
		#QQPlot_Adj(common.rare$results$P.value, common.rare$results$MAP)
		#dev.off()
		
		write.table(assoc$results, file='{output.assoc}', row.names=F, quote=F)
		write.table(common.rare$results, file='{output.cr}', row.names=F, quote=F)
		""")

		
rule generateSSD:
	input	: bed	= BEDFOLDER + '{dataset}.chr{chr}.bed', 
			  bim	= BEDFOLDER + '{dataset}.chr{chr}.bim', 
			  fam	= BEDFOLDER + '{dataset}.chr{chr}.fam',
			  setid	= SSDFOLDER + '{dataset}.chr{chr}.SetId'
	output	: ssd=SSDFOLDER + '{dataset}.chr{chr}.SSD', info=SSDFOLDER + '{dataset}.chr{chr}.Info'
	run		:
		R("""

		library("SKAT")

		bed 	<- '{input.bed}'
		bim 	<- '{input.bim}'
		fam 	<- '{input.fam}'
		setid 	<- '{input.setid}'
		
		ssd 	<- '{output.ssd}'
		info	<- '{output.info}'
		
		Generate_SSD_SetID(bed, bim, fam, setid, ssd, info)
			
		""")

rule generatePhenoSSD:
	input	: bed	= 'pheno_bed/{dataset}.chr{chr}.bed', 
			  bim	= 'pheno_bed/{dataset}.chr{chr}.bim', 
			  fam	= 'pheno_bed/{dataset}.chr{chr}.fam',
			  setid	= PSSDFOLDER + '{dataset}.chr{chr}.SetId'
	output	: ssd=PSSDFOLDER + '{dataset}.chr{chr}.SSD', info=PSSDFOLDER + '{dataset}.chr{chr}.Info'
	run		:
		R("""

		library("SKAT")

		bed 	<- '{input.bed}'
		bim 	<- '{input.bim}'
		fam 	<- '{input.fam}'
		setid 	<- '{input.setid}'
		
		ssd 	<- '{output.ssd}'
		info	<- '{output.info}'
		
		Generate_SSD_SetID(bed, bim, fam, setid, ssd, info)
			
		""")

		
rule generateSetIds:
	input	: set='plinkset/{dataset}.chr{chr}.set'
	output	: setid=SSDFOLDER + '{dataset}.chr{chr}.SetId'	
	script	: 'scripts/generateSetIds.R'

rule generatePhenoSetIds:
	input	: set='pheno_plinkset/{dataset}.chr{chr}.set'
	output	: setid=PSSDFOLDER + '{dataset}.chr{chr}.SetId'	
	script	: 'scripts/generateSetIds.R'
	
	
rule generatePlinkSets:
	input	: expand(BEDFOLDER + '{{dataset}}.chr{{chr}}.{ext}', ext=BED_FILES), geneList='meta/gene.list'	
	output  : 'plinkset/{dataset}.chr{chr}.set', temp('plinkset/{dataset}.chr{chr}.log')
	params 	: bfile=BEDFOLDER + '{dataset}.chr{chr}', out='plinkset/{dataset}.chr{chr}'
	log		: 'logs/{dataset}.chr{chr}.set.log'
	shell	: PLINK 
			+ ' --bfile {params.bfile}'
			+ ' --make-set {input.geneList}'
			+ ' --make-set-collapse-group' 
			+ ' --write-set '
			+ ' --out {params.out} > {log} 2>&1'

	
rule generatePhenoPlinkSets:	
	input	: expand('pheno_bed/{{dataset}}.chr{{chr}}.{ext}', ext=BED_FILES), geneList='meta/gene.list'
	output  : 'pheno_plinkset/{dataset}.chr{chr}.set', temp('pheno_plinkset/{dataset}.chr{chr}.log')
	params 	: bfile='pheno_bed/{dataset}.chr{chr}', out='pheno_plinkset/{dataset}.chr{chr}'
	log		: 'logs/{dataset}.chr{chr}.phenoset.log'
	shell	: PLINK 
			+ ' --bfile {params.bfile}'
			+ ' --make-set {input.geneList}'
			+ ' --make-set-collapse-group' 
			+ ' --write-set '
			+ ' --out {params.out} > {log} 2>&1'
			
rule removeMissingPhenos:
	input	: bed=expand(BEDFOLDER + '{{dataset}}.chr{{chr}}.{ext}', ext=BED_FILES), phenosToRemove=PHENOSTOREMOVE
	output	: bed=expand('pheno_bed/{{dataset}}.chr{{chr}}.{ext}', ext=BED_FILES), tmp='pheno_bed/{dataset}.chr{chr}.log'
	params 	: bfile=BEDFOLDER + '{dataset}.chr{chr}', out='pheno_bed/{dataset}.chr{chr}'
	shell	: PLINK
			+ ' --bfile {params.bfile}'
			+ ' --keep {input.phenosToRemove}'
			+ ' --make-bed'
			+ ' --out {params.out}'

rule generateGeneList:
	input	: refGene=REFGENE
	output	: geneList='meta/gene.list'
	script	: 'scripts/generateGeneList.R'