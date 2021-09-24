# CLoMAT: Conditional Logistic Model Association Tests

## 1. Installation

CLoMAT includes three rare-variant association tests for matched case-control data under the conditional logistic regression (CLR) framework, namely CLR-Burden, CLR-SKAT, and CLR-MiST, as well as a heuristic and fast matching algorithm. CLoMAT provides a general solution to control for population stratification by matching cases and controls based on their ancestry background. It is useful to empower genetic association studies in the setting with a large number of common controls. 

We  built CLoMAT upon the [EPACTS (Efficient and Parallelizable Association Container Toolbox)](https://genome.sph.umich.edu/wiki/EPACTS), which depends on *R* (version 3.2.0) and *groff* (version 1.22.3). After installing EPACTS, we should add the following Rscripts 

- group.clrburden.R
- group.clrskat.R
- group.clrmist.R

to the EPACTS path `./epacts/share/EPACTS/` and install the R packages *survival* and *CompQuadForm* in `./R-3.2.0/lib/`.

The script `heuristic_pairmatch.R`  can be placed in the working directory. The *optmatch* *R* package is required to perform full match. Please install *optmatch* in *R* (version 3.6.0). 

##### Environment variables

- **WK_DIR**: the working directory for running tests.
- **EPACTS_DIR**: the installation directory of EPACTS. 

## 2. Example

### 2.1  Download example data files

For illustration purpose only, we created an example dataset based on [the 1000 Genomes Project]([1000 Genomes | A Deep Catalog of Human Genetic Variation (internationalgenome.org)](https://www.internationalgenome.org/)):

- **EUR.geno.vcf.gz** Genotype file in VCF format for variants extracted from exonic and intronic regions of 503 Europeans. 

- **EUR.pheno.cov.PC.ped** Phenotype and covariate file following the format of  [EPACTS](https://genome.sph.umich.edu/wiki/EPACTS). GROUP indicates case/control status, PC1 and PC2 are genetic ancestry principal components (PCs). In this example, we randomly assigned 100 individuals as cases (GROUP=1) and the rest as controls (GROUP=0). Format of **EUR.pheno.cov.PC.ped**:

  ```
  #FAM_ID	IND_ID	GROUP SEX  PC1	        PC2	  
  HG00114	HG00114	0	  1	   -0.571489	11.9994 
  HG00102	HG00102	1	  2	   2.36762	    21.9064 
  HG00127	HG00127	0	  2	   5.40838	    14.41   
  ```

### 2.2  Perform matching

To control for population stratification, we first match cases and controls based on their ancestry background, such as top PCs. In this example, we matched samples based on PC1 and PC2. We used our heuristic algorithm in `heuristic_pairmatch.R` to perform *1-to-m* match and  the *optmatch* *R* package to perform full match. 

The following commands illustrate how to perform 1-to-1 match:

```R
library(optmatch);
source("${WK_DIR}/Rscripts/heuristic_pairmatch.R");

merged.ped = read.table("${WK_DIR}/data/EUR.pheno.cov.PC.ped",header=T)

match.out = match_on(GROUP ~ PC1+PC2, data = merged.ped,method="euclidean");
match_grp = pairmatch.heuristic(match.out, controls=1, data=merged.ped, width=Inf);
matched.ped = cbind(merged.ped,match_grp=match_grp$matched.group);

IndexNA = which(is.na(match_grp$matched.group));
if(length(IndexNA)>0){matched.ped = matched.ped[-IndexNA,]};

matched.ped = matched.ped[order(matched.ped[,"match_grp"],matched.ped[,"GROUP"]),];

Factor = unique(as.character(matched.ped[,"match_grp"]));
Factor = cbind(Factor,seq(1,length(Factor)));
Index = match(matched.ped[,"match_grp"],Factor[,1]);
matched.ped[,"match_grp"] = Factor[Index,2];

colnames(matched.ped)[1:2]<-c("#FAM_ID","IND_ID");
write.table(matched.ped, file="${WK_DIR}/data/EUR.heuristic.matched.ped",quote=F,row.names=F,sep="\t")
```

The *1-to-m* match is mainly implemented by function 

```
pairmatch.heuristic(D, controls=1, data=NULL, width=Inf)
```

- D: Matching distance matrix. 
- controls: The number *m* of controls to be matched to each case.
- data: The data used for matching. If *NULL*, the output is coded by the row number of matching distance matrix.
- width: Caliper, the width used to exclude case-control pairs with values greater than the width; defaults to *Inf*.

The 1-to-3 match can be done by: 

```R
match_grp = pairmatch.heuristic(match.out, controls=3, data=merged.ped, width=Inf);
```

 The full match can be done by: 

```
match_grp = fullmatch(match.out,data=merged.ped,min.controls=(1/10),max.controls=10);
```

The output file include a new column indicating the matched group: 

```
#FAM_ID	IND_ID	GROUP SEX  PC1	        PC2	      match_grp
HG00114	HG00114	0	  1	   -0.571489	11.9994   1
HG00102	HG00102	1	  2	   2.36762	    21.9064   1
HG00127	HG00127	0	  2	   5.40838	    14.41     2
```

### 2.3 Perform association tests 

To run EPACTS, we need to change the environment path of R to its old versions (R-3.2.0). 

Please refer to the user manual of  [EPACTS](https://genome.sph.umich.edu/wiki/EPACTS) for detailed information of the command lines and output files.

#### 2.3.1 Annotation (EPACTS-3.2.6)

The functional impacts of each variants can be annotated using EPACTS (or other software programs):

```bash
${EPACTS_DIR}/bin/epacts anno \
--in ${WK_DIR}/data/EUR.geno.vcf.gz \
--out ${WK_DIR}/data/EUR.geno.anno.vcf.gz
```

 #### 2.3.2 Define SNP sets (EPACTS-3.2.6)

We define SNP sets using the `make-group` function in EPACTS. The following command generate SNP sets that include missense and nonsense variants in each gene:

```bash
${EPACTS_DIR}/bin/epacts make-group \
--vcf ${WK_DIR}/data/EUR.geno.anno.vcf.gz \
--out ${WK_DIR}/data/EUR.geno.anno.grp \
--format epacts --type Nonsynonymous --type Essential_Splice_Site --type Normal_Splice_Site --type Start_Loss --type Stop_Loss --type Stop_Gain 
```

 #### 2.3.3 Perform gene-based tests (EPACTS-3.2.6; R-3.2.0)

In the following examples, we perform gene-based tests. We filter variants with  minor allele frequency MAF<0.01, and adjust for sex as a covariate. CLR-based test can be specified by `--test`. The weights of CLR-based tests can be specified by `--beta beta1,beta2`, which means the weights follow the beta distribution with Beta(MAF, beta1, beta2). The strata information from matching is specified as the *last* covariate (`--cov`). The number of parallel jobs is set by `--run`. More options can be found on the website of  [EPACTS](https://genome.sph.umich.edu/wiki/EPACTS).

To run **CLR-Burden**:

```bash
${EPACTS_DIR}/bin/epacts group \
--vcf ${WK_DIR}/data/EUR.geno.anno.vcf.gz --group ${WK_DIR}/data/EUR.geno.anno.grp \
--ped ${WK_DIR}/data/EUR.heuristic.matched.ped \
--out ${WK_DIR}/output/test.clrburden \
--test clrburden --beta 1,25 \
–-min-maf 0.01 –-max-maf 0.5 \
--pheno GROUP --cov SEX --cov match_grp \
--run 2
```

The output includes 11 files of individual information, covariates information and results. 

The manhattan and QQ plots are in `output/*.epacts.mh.pdf` and `output/*.epacts.qq.pdf`. 

The test results are stored in `output/*.epacts`, and the top 5000 association signals are output to `output/*.epacts.top5000`. Below illustrates the example output of  **CLR-Burden**. The first 9 columns follow the standard output of [EPACTS](https://genome.sph.umich.edu/wiki/EPACTS).

```
#CHROM BEGIN END MARKER_ID NS FRAC_WITH_RARE NUM_ALL_VARS NUM_PASS_VARS	NUM_SING_VARS PVALUE LRT_STAT SCR_PVAL SCR_STAT
 1  77752622  78024332      1:77752622-78024332_AK5 200 0.435 28  2  1 0.82243 0.05036700 0.82309 0.04998300
 1 109358807 109506029 1:109358807-109506029_AKNAD1 200 0.675 76  8  6 0.46843 0.52568000 0.46863 0.52522000
```

- NS: Number of phenotyped samples with non-missing genotypes.
- FRAC_WITH_RARE: Fraction of individuals carrying variants  passing the `--min-maf`, `--max-maf`, `--min-mac` thresholds.
- NUM_ALL_VARS: Number of all variants defining the group.
- NUM_PASS_VARS: Number of variants passing the `--min-maf`, `--max-maf`, `--min-mac` thresholds.
- NUM_SING_VARS: Number of singletons among all variants.
- PVALUE: *P* value from likelihood ratio test.
- LRT_STAT: Q statistic from likelihood ratio test.
- SCR_PVAL: *P* value from score test.
- SCR_STAT: Q statistic from score test.

Similarly, **CLR-SKAT** and **CLR-MiST** can be performed by changing `--test clrburden`  to `--test clrskat` and `--test clrmist` in above command, respectively.

The example output of **CLR-SKAT**:

```
#CHROM BEGIN END MARKER_ID NS FRAC_WITH_RARE NUM_ALL_VARS NUM_PASS_VARS	NUM_SING_VARS PVALUE LRT_STAT SCR_PVAL SCR_STAT
 1  77752622  78024332      1:77752622-78024332_AK5 200 0.435 28  2  1 0.81658  0.23294 0.81658  0.23294
 1 109358807 109506029 1:109358807-109506029_AKNAD1 200 0.675 76  8  6 0.56411 13.64900 0.56751 13.56700
```

The example output of **CLR-MiST**:

```
#CHROM BEGIN END MARKER_ID NS FRAC_WITH_RARE NUM_ALL_VARS NUM_PASS_VARS	NUM_SING_VARS PVALUE CHISQ_STAT Burden_PVAL	Burden_STAT	SKAT_PVAL	SKAT_STAT
 1  77752622  78024332      1:77752622-78024332_AK5 200 0.435 28  2  1 0.98150 0.41177 0.82243 0.050367 0.98966 1.3879e-07 
 1 109358807 109506029 1:109358807-109506029_AKNAD1 200 0.675 76  8  6 0.53506 3.13770 0.46843 0.525680 0.44466 1.2636e+01 
```

- PVALUE: *P* value from hybrid mixed effects score test.
- CHISQ_STAT: Chi-square statistic in Fisher's procedure to combine CLR-Burden and CLR-SKAT-type test statistics. 
- Burden_PVAL: *P* value from CLR-Burden test.
- Burden_STAT: Q statistic from CLR-Burden test.
- SKAT_PVAL: *P* value from CLR-SKAT test (regressing out genetic burden).
- SKAT_STAT: Q statistic from CLR-SKAT test (regressing out genetic burden).
