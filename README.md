## AIPS (Absolute Inference of Patient Signatures)

This github repository contains the R package that could be used to partition breast cancer gene expression profile(s) using the 1733 models defined in AIPS.

Here is how you could install AIPS in R :

* Download the AIPS package from github (https://github.com/meoyo/AIPS)

* Install AIPS in R using eg. this command from the shell :

`> R CMD INSTALL AIPS_1.0.0.tar.gz`

* The easisest way to know how to run AIPS is to read the vignette by launching R and then load AIPS using :

```R
library(AIPS)
library(Biobase)
openVignette()
```

## AIPS applied on RNA-seq using Kallisto

If you are interested to apply AIPS on RNA-seq samples the easiest way is to follow the instruction in:

kallisto.AIPS/INSTALL

If you have any questions please contact Eric Paquet at (eric.r.paquet@gmail.com)
