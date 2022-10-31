# mev-lrt

This repository contains a WebMeV-compatible tool for performing differential expression analysis using the likelihood ratio test Bioconductor's DESeq2. 

Note that this assumes you are testing the significance of a single covariate, such as an applied treatment. In R's model syntax, we are testing `~treatment` against a null model of `~1` (not controlling for any other covariates). This could also be used, for instance, to test whether `batch` is a significant factor in a given experiment. Note, however, that we do not check other covariates for confounding.

The outputs are:
- A tab-delimited file of the differential expression results merged with the normalized counts
- A tab-delimited file of just the normalized counts.
- A string telling us which inter-group comparison produced the log2 fold change estimate.
- A json-format string giving the mapping of groups to samples. For example,
```
{
    "groupA": ["sampleA", "sampleB"],
    "groupB": ["sample1", "sample2"]
}
```
(This is mainly for convenience for our WebMeV frontend. The annotations come directly from the annotation file provided as an input.)

The concatenation of the differential expression results and counts is for convenience with the WebMeV frontend interface since it avoids pulling data from two different files.

Since the likelihood ratio test is capable of comparing $n>=2$ groups, the tests (and resulting p-values) correspond to testing the significance of the covariate of interest. For fold-change estimates, there are ${n}\choose{2}$ choices, and DESeq2 will report one of those. We provide that as one of the outputs so that it's clear what the fold-change estimate corresponds to.

---

### To run external of WebMeV:

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/deseq2:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-lrt/pkgs/container/mev-lrt)

To run, enter the container in an interactive shell:
```
docker run -it -v$PWD:/work <IMAGE>
```
(here, we mount the current directory to `/work` inside the container)


Then, run the script:
```
Rscript /opt/software/deseq2.R \
    <path to raw/integer counts> \
    <path to annotations> \
    <covariate/column name>
```
The call to the script assumes the following:
- The input file of expression counts is tab-delimited format and contains only integer entries
- The annotation matrix contains the sample IDs in the first column (such that they have *some* intersection with the sample names in the count matrix) and other covariates in the following columns (e.g. phenotype, treatment, experimental group, etc.)
- The final argument providing the column name *exactly* matches that column header in the annotation file.