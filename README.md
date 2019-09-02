# GWAS-pipeline

## Description

This is a suit of GWAS workflows for [CHIMGEN]([http://chimgen.tmu.edu.cn](http://chimgen.tmu.edu.cn/)) project built on [NEXTFLOW](<https://www.nextflow.io/>) framework.

## Contents

- [Getting started](#4)

- [QC of raw data](#1)

- [Imputation](#2)

- [Association analysis](#3)

## <a name="4"></a>Getting started

I built the whole pipeline on [NEXTFLOW](<https://www.nextflow.io/>) framework. And NEXTFLOW can be installed by conda. To get started, you should clone this repo firstly.

```shell
git clone https://github.com/Jianhua-Wang/GWAS-pipeline.git
```

then set up the conda environment:

```shell
conda env create -f ./bin/gwascondaenv.yml
```

and use the commands below to activate/deactivate the environment:

```shell
conda activate gwas
conda deactivate
```

This usage of this pipeline is quite simple, for example:

```shell
nextflow run qc.nf -c qc.config
```

`qc.nf` is the NEXTFLOW script for GWAS quality control and `qc.config` is the configuration file that stores the parameters of the pipeline. NEXTFLOW will genere a `work` folder in work directory for debuging and storing intermediate files. You can remove the Intermediate files by manually or using `nextflow clean -f`.

---

## <a name="1"></a>QC of raw data

### Purpose

This is a pipeline of genotype data quality assessment and control, including variant-level QC and individual-level QC. These steps are critical for the following imputation and association analysis, could reduce the false positive of GWAS results.

I referred to the *Nature Protocols* from [Anderson et al.](<https://www.nature.com/articles/nprot.2010.116>) and a NEXTFLOW workflow from [H3ABioNet](<https://github.com/h3abionet/h3agwas>). I didn't describe the details of every step in this tutorial, but it's very easy to understand the principles on the *Nature Protocols*.

### Input

PLINK binary .bed, .bim, and .fam file

There are various types of genotype formats, such as PLINK bed, gen, vcf, bgen, even pgen. You can find the description of them on [PLINK's website](https://www.cog-genomics.org/plink/1.9/formats). In this part, we mainly use PLINK bed format, which contains a [bed](https://www.cog-genomics.org/plink/1.9/formats#bed) file for genotype, a [bim](https://www.cog-genomics.org/plink/1.9/formats#bim) file for variant information, and a [fam](https://www.cog-genomics.org/plink/1.9/formats#fam) file for individual information, with the same prefix.

### Output

This pipeline does not only output the clean data but also generate tables of filtered variants and individuals and figures of data assessment. I put the tables and figures in a HTML report to read easily.

| Output file              | Description                                       |
| ------------------------ | ------------------------------------------------- |
| test-nd-c-c.bim          | bim file of clean data                            |
| test-nd-c-c.bed          | bed file of clean data                            |
| test-nd-c-c.fam          | fam file of clean data                            |
| test.dups                | duplicated variants                               |
| test.badsex              | sample with ambiguous sex                         |
| test-nd-c.irem           | removed sample due to missing rate                |
| test-nd-fail_IBD.txt     | related sample                                    |
| test-nd-c-fail_het.txt   | sample with extreme high or low heterozygous rate |
| test-nd-c-c.irem         | sample removed in phased 3                        |
| test-GWAS-QC_report.html | report of QC                                      |
| test_runtime.html        | report of workflow running                        |
| test_timeline.html       | report of workflow running                        |

### Usage

