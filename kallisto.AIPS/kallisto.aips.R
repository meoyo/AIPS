## Usage :
## Rscript kallisto.aips.R paired.end.1.fastq.gz paired.end.2.fastq.gz output_folder
##
## Output:
##
## All the standard kallisto output in output_folder
##
## AND the AIMS and AIPS partitions in:
##
## output_folder/AIMS.AIPS.xls
##
## --------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------
## INSTALL
## 
## To be able to run this script you need to install those
## R/Bioconductor packages :
## 
## AIMS : https://www.bioconductor.org/packages/release/bioc/html/AIMS.html
## AIPS : https://github.com/meoyo/AIPS
## org.Hs.eg.db : https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
##
## You also need to install kallisto :
##
## https://pachterlab.github.io/kallisto/starting.html
##
## and a transcriptomes:
##
## http://bio.math.berkeley.edu/kallisto/transcriptomes/
##
## this example will work with the Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz file.
##
## You might have to update the variable KALLISTO.COMMAND to fit your configuration.
##
## If you have any questions or concerns contact Eric Paquet at (eric.r.paquet@gmail.com)
##
## Copyright McGill University 2016
##
start.time=proc.time()
KALLISTO.COMMAND = "kallisto quant -i Homo_sapiens.GRCh38.rel79.cdna.all.idx --plaintext -o %s %s %s"

args <- commandArgs(trailingOnly = T)

stopifnot(length(args) == 3)

paired.end.1 <- args[1]
paired.end.2 <- args[2]
output <- args[3]

## output <- sprintf("kallisto.aips.%d.out",Sys.getpid())

message(sprintf("## RUNNING KALLISTO AIPS ON :\n## Paired-end file 1 = %s\n## Paired-end file 2 = %s\n## output in = %s\n\n",
                paired.end.1,paired.end.2,output))

system(sprintf(KALLISTO.COMMAND,output,paired.end.1,paired.end.2))

require(AIPS)
require(AIMS)
require(org.Hs.eg.db)

options(stringsAsFactors=F)
kallisto.tpm <- read.delim(file.path(output,"abundance.tsv"),header=T,sep="\t")
ensembl2entrez <- select(org.Hs.eg.db, kallisto.tpm[,1], "ENTREZID", "ENSEMBLTRANS")
ensembl2entrez <- ensembl2entrez[match(kallisto.tpm[,1],ensembl2entrez[,1]),]
aims <- apply.AIMS(matrix(kallisto.tpm[,"tpm"],ncol=1),as.character(ensembl2entrez[,2]))
aips <- apply.AIPS(matrix(kallisto.tpm[,"tpm"],ncol=1),as.character(ensembl2entrez[,2]))

aips.out <- cbind(aips$gs.info,AIPS.ASSIGNMENT=aips$cl,AIPS.POSTERIOR=aips$posterior)
elapse.time <- (proc.time() - start.time)[3]

to.w <- rbind(c("Paired-end 1",paired.end.1,"elapse-time (seconds)",elapse.time,""),
              c("Paired-end 2",paired.end.2,"","",""),
              c("AIMS-subtype","Paquet et al. JNCI","",aims$cl,aims$posterior),
              aips.out)

colnames(to.w)[4:5] <- c("Assignments","Posterior")

out.file <- file.path(output,"AIMS.AIPS.xls")
message(sprintf("## Writing AIMS and AIPS assignments in %s",out.file))

write.table(to.w,sep="\t",col.names=T,row.names=F,file=out.file)
