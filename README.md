# XWAS Cleaning Pipeline

In this GitHub we implemented almost all steps from the article "Common X-Chromosome Variants Are Associated with Parkinson Disease Risk" (https://doi.org/10.1002/ana.26051)

We did not implemented the steps based on ancestry and the steps related to the imputation panel (which are performed by TopMed).

The pipeline has cleaning process to Autosomal and X-chromosome.

## Autosomal

The cleaning process to autosomal is composed by:

- Remove individuals without covariatives
  - Mandatory covariates: Sex and Phenotypes
  - All other coavariates in the covar file will be checked and those samples that have "NA" will be removed
- Remove SNPs located in Structural variants (file TryTyper.txt made by Valerio Napolioni)
- Remove duplicated SNPs
- Remove monomorphic SNPs
- Remove potential probe sites using gnomAD. 

To run this step we downloaded the gnomAD data (gnomad.genomes.v3.1.1.sites.chr*.vcf.bgz) and we kept just the columns with AF and AF > 1%
```
bcftools annotate -x ^INFO/AF -Oz -o gnomadOnlyAF_1.vcf.gz gnomad.genomes.v3.1.1.sites.chr1.vcf.bgz --threads 10 
bcftools index gnomadOnlyAF_1.vcf.gz --threads 10
bcftools view -i "AF > 0.01" gnomadOnlyAF_1.vcf.gz -Oz -o gnomadOnlyAF_onePercent_1.vcf.gz
bcftools index gnomadOnlyAF_onePercent_1.vcf.gz --threads 10
```
- Remove SNPs with HWE in controls with p-value lower than 10^-5
- Remove SNPs with missing data in SNPs grater than 5%
- Relationship control using king to calculate and NAToRA to remove (https://github.com/ldgh/NAToRA_Public)

## Chromosome X

The cleaning process to autosomal is composed by:

- Keep the remaining individuals on Autosomal cleaning process
- Remove SNPs located in Structural variants (file TryTyper.txt made by Valerio Napolioni)
- Remove duplicated SNPs
- Remove monomorphic SNPs
- Remove potential probe sites using gnomAD. 
- Remove SNPs with HWE in controls lower than 10^-5
- Remove SNPs with missing data in SNPs grater than 5%
- Remove individuals that fail in sex-check
- Heterozygous SNPs was converted to missing data in males
- Remove SNPs with differential missingness between cases and controls with p-values lower than 10^-5
- Remove SNPs with differential missingness between males and females with p-values lower than 10^-5
- Remove SNPs with differential MAF between males and females with p-values lower than 10^-5
- Phasing using Eagle v2.4.1 with 1KGP 30x as reference to swap REF/ALT alleles (--allowRefAltSwap)

## Programs used

In this pipeline we use:

- Plink 
- Python3
  - NetworkX (NAToRA) *
  - Scipy (Fisher's exact test)
  - Numpy
- King *
- NAToRA *

- BGZIP **
- BCFTools **
- Eagle **

You can change the path of the programs using the respective flags for each software

\* Used on relationship control

\*\* Used swap reference and alternative alleles using 1KGP 30x as reference
## Options 

```
usage: harmonizationAdapted.py [-h] -a AUTOSOMAL -x XCHROMOSOME -s STRUCTURAL -o OUTPUT -O ONETHOUSAND -f FOLDER -c COVAR -g GENETICMAP -p PHENOTYPE [-I INTERPRETER] [-P PLINK] [-T PLINK2] [-B BCFTOOLS] [-b BGZIP] [-E EAGLE]
                               [-N NATORA] [-K KING] [-d DEGREE] [-C] [-G GNOMAD] [-r REFERENCEBUILD] [-t THREADS] [-H]

XWAS cleaning

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -a AUTOSOMAL, --autosomal AUTOSOMAL
                        Autosomal input file
  -x XCHROMOSOME, --xchromosome XCHROMOSOME
                        Chromosome X input file
  -s STRUCTURAL, --structural STRUCTURAL
                        Variants in strutural region
  -o OUTPUT, --output OUTPUT
                        Output file prefix
  -O ONETHOUSAND, --oneThousand ONETHOUSAND
                        Path for One Thousand Genomes vcf file
  -f FOLDER, --folder FOLDER
                        Folder to store output files
  -c COVAR, --covar COVAR
                        File with the covariatives
  -g GENETICMAP, --geneticMap GENETICMAP
                        Path to genetic Map from eagle
  -p PHENOTYPE, --phenotype PHENOTYPE
                        Phenotype column on covar file

Optional arguments:
  -I INTERPRETER, --interpreter INTERPRETER
                        Path to python3 (default: python)
  -P PLINK, --plink PLINK
                        Plink path (default: plink)
  -T PLINK2, --plink2 PLINK2
                        Plink2 path (default: plink2)
  -B BCFTOOLS, --bcftools BCFTOOLS
                        BCFTools path (default: bcftools)
  -b BGZIP, --bgzip BGZIP
                        Bgzip path (default: bgzip)
  -E EAGLE, --eagle EAGLE
                        Eagle path (default: eagle)
  -N NATORA, --NAToRA NATORA
                        NAToRA path (default: NAToRA)
  -K KING, --king KING  KING path (defaul: king)
  -d DEGREE, --degree DEGREE
                        Relationship degree to be removed. (default = "", that meansto not remove relationship)
  -C, --continuous      Set the phenotype as continuos
  -G GNOMAD, --gnomAD GNOMAD
                        Path to gnomAD files with the chr replaced by * (default = "", not perform this step)
  -r REFERENCEBUILD, --referenceBuild REFERENCEBUILD
                        Human Genome build (default : hg38)
  -t THREADS, --threads THREADS
                        Number of threads to Eagle (default = 1)
  -H, --homozygousMale  The males has homozygous diploid chr X (default False)
```

## Example of command line

```
python3.8 harmonizationAdapted.py \
-a /home/user/input/Input_Autosomal \
-x /home/user/input/Input_chrX  \
-s /home/user/input/TriTyper.txt \
-o OutputTest \
-c /home/user/input/Inputcovar.txt \
-p DISEASE \
-f /home/user/input/OutTest \
-G /home/user/input/gnomAD/gnomadOnlyAF_onePercent_\*.vcf.gz \
-d 2 \
-N ./NAToRA_Public.py \
-I python3.8 \
-t 40 \
-K ./king \
-g /home/user/input/genetic_map_hg38_withX.txt.gz \
-O /home/user/input/OneThousand/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased_header.vcf.gz
```

# XWAS PCA and running

After the harmonization process, we imputed our data using TOPMED imputation panel. After this we perform the regression using the script makeAllStepsToRegression.py

```
usage: makeAllStepsToRegressionV2.py [-h] [-g GENOTYPED] [-i IMPUTED] -t TABLECOVAR -C COUNTRYFILE [-c COUNTRY [COUNTRY ...]] -l COVARLIST [COVARLIST ...] -m MAXPC -r R2 [-F] [-H] [-n NAME] [-f FOLDER] [-G GWAMA] [-p PLINK2]
                                     [-P PYTHON]

PCA and regression

optional arguments:
  -h, --help            show this help message and exit

Data arguments:
  -g GENOTYPED, --genotyped GENOTYPED
                        Genotyped file name
  -i IMPUTED, --imputed IMPUTED
                        Imputed file name
  -t TABLECOVAR, --tableCovar TABLECOVAR
                        File with covariatives to be added to the model
  -C COUNTRYFILE, --countryFile COUNTRYFILE
                        File with relation Ind country
  -c COUNTRY [COUNTRY ...], --country COUNTRY [COUNTRY ...]
                        Country to analyze (default: all). You can select more than one country

Parameter arguments:
  -l COVARLIST [COVARLIST ...], --covarList COVARLIST [COVARLIST ...]
                        List of covar to be used (do not include PCAs)
  -m MAXPC, --maxPC MAXPC
                        max PC to be used on covar
  -r R2, --r2 R2        r2 cutoff
  -F, --firth           Force all PLINK2 regressions use the firth
  -H, --homozygousOnly  Remove heterozygous from imputed file

Output arguments:
  -n NAME, --name NAME  Name to use
  -f FOLDER, --folder FOLDER
                        Folder to output files

Programs arguments:
  -G GWAMA, --gwama GWAMA
                        GWAMA program (default = gwama)
  -p PLINK2, --plink2 PLINK2
                        Path of PLINK 2 (default = plink2)
  -P PYTHON, --python PYTHON
                        Path of Python 3 (default = python)

```

### Example of command line

```
python3.8 makeAllStepsToRegression.py \
-g /home/user/input/OutTest/OutputTest_chrX_Phased.vcf.gz \
-i /home/user/input/OutTest/OutputTest.chrX.imputed.dose.vcf.gz \
-n OutputTest \
-f /home/user/input/OutTest/ResultOutTest \
-t /home/user/input/Inputcovar.txtt \
-C /home/user/input/IDCountry.txt \
-l AGE SEX \
-m 10 \
-r 0.8 -F -H \
-p ./plink2 \
-P python3.8 \
-G ./GWAMA \
```
The parameter -F forces [PLINK2 use the firth regression](https://www.cog-genomics.org/plink/2.0/assoc#no_firth), that is good for small sample size. 

The parameter -H was a test which we make all heterozygous variants for females being missing data. The idea behind this test is based on the fact that one of the two X chromosomes is silenced in females, so there is no way to say *in silico* which of the chromosomes is used for biological processes. So, the heterozygous could act like a noise to regression, since there is no way to predict which allele was used on transcription. In our tests, this flag gave more statistical power to females analysis, but caused a loss of power to analysis which males and females are analyzed together.

We also uploaded the "makeAllStepsToRegressionAutosomalPCA.py", that perform the same analysis than "makeAllStepsToRegression.py" but using the Autosomal PCA. This changed was made based on the reviewer suggestion. This script will ask the autosomal dataset (-a) and we recomend to use the "<name>_Relatedness" from the harmonization pipeline.

 ## Acknowledgements

This work is supported by NIH Grant R01 1R01NS112499-01A1, MJFF Grant ID: 18298, ASAP-GP2 and Parkinson's Foundation

 The authors want to thank Dr Marla Mendes de Aquino (Genetics and Genome Biology Program, The Hospital for Sick Children, Toronto, ON, Canada) for the suggestions to improve the pipeline
 ### Contact
 
 Developer: Thiago Peixoto Leal. PhD (PEIXOTT@ccf.org or thpeixotol@hotmail.com)


