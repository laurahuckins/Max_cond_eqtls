# Max_cond_eqtls

Files in this repository can be used to iteratively calculate max and conditional eqtls.

I have iterated over 10 differents sets if IDs as well, in order to be compatible with a later 10-fold cross validation step. To use all ids, simply remove those lines from the code. (Commented with :remove for all ids:)

This folder contains the following code snippets:

1. Get_files.R :

This file is a messy way of getting all the files I need in the right format, and removing any unwanted ids from the list. I take no responsibility for this file!
The outputs of this file are:

vcf_file: this contains the dosages for all snps that you want to test. Rownames= snps, colnames= samples ids
rnaseqfile: this contains gene expr values for all genes. Rownames= genes, colnames= samples ids, same order as vcf_file
meta: my covariate file. 
snpspos: a data.frame file, colnames(snpspos)<-c("snpid", "chr", "pos")
genepos: a data.frame file, colnames(genepos)<-c("geneid", "chr", "left", "right")

2. get_max_eqtls.R

This file takes the output of Get_files.R and creates a formatted maxeqtlt file. If you only want maxeqtls, you're done now! This will give you a max eqtlt for every gene. 


