# This file takes the output of Get_files.R and creates a formatted maxeqtlt file. If you only want maxeqtls, you're done now!


useModel = modelLINEAR

cvrt = SlicedData$new();
cvrt$CreateFromMatrix(as.matrix(t(cov)))

gene = SlicedData$new();
gene$CreateFromMatrix(as.matrix(rnaseqfile))
rownames(gene)<-rownames(rnaseqfile)

snps=SlicedData$new()
snps$CreateFromMatrix(as.matrix(vcf_file))

useModel=modelLINEAR
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
#output_file_name     = file.path(outdir, transfilename),
pvOutputThreshold     = 0,
useModel = useModel,
#errorCovariance = errorCovariance,
verbose = TRUE,
#output_file_name.cis = file.path(outdir, cisfilename),
output_file_name.cis = "temp.cis",
pvOutputThreshold.cis = 0.1,
snpspos = snpspos,
genepos = genepos,
cisDist = 1e6,
#pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

# Now: get list of max eqtls:

eqtls=me$cis$eqtls

# get max snp per gene

min_pvals=me$cis$min.pv.gene

genes=unique(eqtls$gene)

assign_maxeqtl<-function(gene, eqtls){

maxeqtl<-eqtls[grep(gene,eqtls$gene)[1],]

return(maxeqtl)

}

maxeqtllist.tmp<-lapply(genes, assign_maxeqtl, eqtls)

maxeqtllist<-do.call(rbind, maxeqtllist.tmp)

rownames(maxeqtllist)<-maxeqtllist$gene
