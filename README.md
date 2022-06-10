# XWAS Cleaning Pipeline

In this GitHub we implemented almost all steps from the article "Common X-Chromosome Variants Are Associated with Parkinson Disease Risk" (https://doi.org/10.1002/ana.26051)

We did not implemented the steps based on ancestry and the steps related to the imputation panel (which are performed by TopMed).

The pipeline has cleaning process to Autosomal and X-chromosome.

## Autosomal

The cleaning process to autosomal is composed by:

- Remove individuals without covariatives (Age, Sex and Phenotype)
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
usage: main.py [-h] -a AUTOSOMAL -x XCHROMOSOME -s STRUCTURAL -o OUTPUT -O ONETHOUSAND -f FOLDER -c COVAR -g GENETICMAP -p PHENOTYPE [-I INTERPRETER] [-P PLINK] [-T PLINK2] [-B BCFTOOLS] [-b BGZIP] [-E EAGLE] [-N NATORA] [-K KING]
               [-d DEGREE] [-C] [-G GNOMAD] [-r REFERENCEBUILD] [-t THREADS]

XWAShilpa

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
                        Relationship degree to be removed. (default = "", that means to not remove relationship)
  -C, --continuous      Set the phenotype as continuos
  -G GNOMAD, --gnomAD GNOMAD
                        Path to gnomAD files with the chr replaced by * (default = "", not perform this step)
  -r REFERENCEBUILD, --referenceBuild REFERENCEBUILD
                        Human Genome build (default : hg38)
  -t THREADS, --threads THREADS
                        Number of threads to Eagle (default = 1)
```

## Example of command line

```
python3.8 main.py \
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




