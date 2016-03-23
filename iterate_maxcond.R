# This file will iteratively find conditional eqtls, conidiotning on max, then max + cond1, etc.

# Inputs: maxeqtl file from get_max.R, initial files loaded from Get_files.R

# This function takes the cross-validation fold you want to test.
# So set=1 means first cross-validation iteration, so the 40 test ids for that iteration will be excluded.

# input.row is the row relating to the gene of interest, from get_max_eqtls.R

iterate<-function(set, genename){

# Get the max eqtl for that gene, for this cv set.

gene.locs<-which(maxeqtllist.sets[[set]][,2]==genename)
max.use<-gene.locs[which(maxeqtllist.sets[[set]][gene.locs,4]==min(maxeqtllist.sets[[set]][gene.locs,4]))]

if(length(max_use)==1){ # If we only have one max eqtl...

# The first section sets us up to deal with set-specific ids. Comment this section to use all, and set rnaseqfile.set=rnaseqfile, etc.

setids<-sets[[set]][,1]

input.row=maxeqtllist.sets[[set]][max.use,]

vcf_file.set<-vcf_file[,-(which(colnames(vcf_file)%in%setids))]
cov.set<-cov[,-(which(colnames(vcf_file)%in%setids))]
rnaseqfile.set<-rnaseqfile[,-(which(colnames(rnaseqfile)%in%setids))]

# End of section to deal with individual sets

iteration =2

genename=as.character(input.row[1,2])
snploc=which(rownames(vcf_file.set)==as.character(input.row[1,1]))

output.row<-cbind(input.row, "1")
colnames(output.row)[7]<-"iteration"

gene1=SlicedData$new();
gene1$CreateFromMatrix(as.matrix(rnaseqfile.set[genename,]))

pass_flag=1


while(pass_flag==1){

cov1=rbind(cov.set, vcf_file.set[snploc,])
cvrt1 = SlicedData$new();
cvrt1$CreateFromMatrix(as.matrix(cov1))

snps1=SlicedData$new()
snps1$CreateFromMatrix(as.matrix(vcf_file.set[-(snploc),]))


useModel=modelLINEAR
me1 = Matrix_eQTL_main(
snps = snps1,
gene = gene1,
cvrt = cvrt1,
#output_file_name     = file.path(outdir, transfilename),
pvOutputThreshold     = 0,
useModel = useModel,
#errorCovariance = errorCovariance,
verbose = TRUE,
#output_file_name.cis = file.path(outdir, cisfilename),
output_file_name.cis = "temp.cis",
pvOutputThreshold.cis = threshold*2,
snpspos = snpspos,
genepos = genepos,
cisDist = 1e6,
#pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

if(me1$cis$eqtls[1,5]<=0.2){
snploc=c(snploc,which(rownames(vcf_file)==as.character(me1$cis$eqtls[1,1])))
new<-cbind(me1$cis$eqtls[1,], as.character(iteration))
colnames(new)<-colnames(output.row)
output.row<-rbind(output.row, new)
} else {
new<-cbind(me1$cis$eqtls[1,], "fail")
colnames(new)<-colnames(output.row)
output.row<-rbind(output.row, new)
#output.row
pass_flag=0 # we no longer reach significance
}

iteration=iteration+1
}
} else {    # what to do if we have two max eqtls

output.row<-maxeqtllist.sets[[set]][max.use,]
output.row<-cbind(output.row,"Mult_maxes")

}

output.row<-cbind(output.row, output.row[,7])
output.row[,8]<-as.character(set)
colnames(output.row)[8]<-"Set"
return(output.row)

}

# The function below iterates over all sets and produces an output file with eqtls listed by cross-validation fold, and max/cond status.

get_nice_output<-function(index, genes){  # This will allow us to iterate over all the sets and get a single output file

# here, iterate over set lists as well
output=mclapply(1:10,iterate, as.character(genes[index]))

# mclappaly is a parrallelised version of lapply; for local environments etc switch to lapply

output.sets=do.call(rbind, out)
output.sets[,1]<-output.sets[1,1] # Make sure all first row info is the genename as character
return(output.sets)

}

# Finally, call both functions in lapply to get a file with iterations over all cross-validation folds, and all max/cond eqtls:

output.new=mclapply(1:3, get_nice_output, genes) # mclapply is // lappaly, see above comment.
output=do.call(rbind, output.new)

