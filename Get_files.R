# This file is a messy way of getting all the files I need in the right format, and removing any unwanted ids from the list. I take no resposnisbility for this file!
# The outputs of this file are:

# vcf_file: this contains the dosages for all snps that you want to test. Rownames= snps, colnames= samples ids
# rnaseqfile: this contains gene expr values for all genes. Rownames= genes, colnames= samples ids, same order as vcf_file
# meta: my covariate file. 
# snpspos: a data.frame file, colnames(snpspos)<-c("snpid", "chr", "pos")
# genepos: a data.frame file, colnames(genepos)<-c("geneid", "chr", "left", "right")



chr=10

library("MatrixEQTL")

# Load in necessary files:

rnaseqfile<-read.delim("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/DLPFC.ensembl.KNOWN_AND_SVA.ADJUSTED.VOOM_NORMALIZED.GE.WEIGHTED_RESIDUALS.chr10.tsv", header=T, sep='')

#rnaseqfile<-read.delim("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/Chr10_exons", header=T, sep='')

vcf_file<-read.delim(sprintf("/sc/orga/projects/CommonMind/scratch/adobbyn/eqtl/CMC_from_SYNAPSE/CM5-chr%s_imputed.dos.id",chr), header=T, sep='')

vcf_file.perm<-vcf_file # keep a copy just in case as this is slow to load

snplist=read.delim("/sc/orga/work/huckil01/filt_file/Chr10.bim", header=F)

# for each gene, cycle through all snps in the chr (that fall within snplist)


dupl<-which(duplicated(vcf_file$SNP)==TRUE)
vcf_file<-vcf_file[-(dupl),]

keep<-which(vcf_file$SNP %in% snplist[,2])

vcf_file<-vcf_file[keep,]
rownames(vcf_file)<-vcf_file$SNP

snpspos<-vcf_file[keep,c(2,1,3)]
snpspos[,2]<-10
snpspos<-data.frame(snpspos)
colnames(snpspos)<-c("snpid", "chr", "pos")

rownames(rnaseqfile)<-rnaseqfile$Gene

sets=vector(10, mode="list")    # This is so I can calc for each set of indivs at the same time

for (i in 1:10){

sets[[i]]<-read.delim(sprintf("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/BSLMM/Ids.set%s", i), header=F, sep='')

}


gencode<-read.delim("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/genfile_eQTLS_cond", header=F, sep='')
rownames(gencode)<-gencode[,9]

genes.use<-intersect(rownames(gencode), rownames(rnaseqfile))

rnaseqfile<-rnaseqfile[genes.use,]
gencode<-gencode[genes.use,]

gencode<-gencode[,1:9]
gencode[,5]<-gencode[,4]-gencode[,3]

genepos<-gencode[,c(9,1,3,4)]
genepos[,2]<-10
genepos<-data.frame(genepos)
colnames(genepos)<-c("geneid", "chr", "left", "right")

rmfile=read.delim("/sc/orga/projects/CommonMind/lhuckins/non-caucasian.ids", header=F, sep='')

vcf_rm<-which(colnames(vcf_file)%in%rmfile[,1])

vcf_file<-vcf_file[,-(vcf_rm)]

indivs.use<-intersect(colnames(rnaseqfile), colnames(vcf_file))

vcf_file<-vcf_file[,indivs.use]
rnaseqfile<-rnaseqfile[,indivs.use]

vcf_file<-as.matrix(vcf_file)
rnaseqfile<-as.matrix(rnaseqfile)

#gene1=SlicedData$new(matrix(rnaseqfile, nrow=1))
#snps1=SlicedData$new(matrix(vcf_file, nrow=1))

#gene1=SlicedData(matrix(rnaseqfile[1,]))
#snps1=SlicedData(matrix(vcf_file[1,]))

meta<-read.delim("/sc/orga/projects/CommonMind/data/MERGED_METADATA.csv", header=T, sep=",")

metaids<-read.delim("/sc/orga/projects/CommonMind/lhuckins/Meta.ids.matched", header=F, sep='')

metaids.up<-toupper(metaids[,1])

meta<-meta[which(metaids.up %in% indivs.use),]
meta.use<-metaids.up[which(metaids.up %in% indivs.use)]

rownames(meta)<-meta.use
