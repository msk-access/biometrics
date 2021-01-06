---
description: Step for preparing the BAM file(s)
---

## Extract
Running this step is a **prerequisite** before running any of the other tools. This step extracts the pileup and coverage information from your BAM file(s) and stores the result in a Python pickle file (which contains JSON data). You can determine where to store the output files by specifying `-db` argument. This allows for much faster analyses that make repeated use of your samples.

There are two main types of required input:

- **Sample information:** the BAM file and any associate annotation (e.g. sample grouping).
- **Supporting files:** reference FASTA, VCF, and BED file.

Moreover, there are two ways to provide the sample information: (1) provide a CSV file, or (2) specify via the command line arguments.

#### CSV input
This method is easier for when you have many samples. Just provide a CSV file with five columns: sample name, sample group, sample type, sample sex, and path to the sample's BAM file. An example with three samples is below:

```text
sample_name,sample_group,sample_type,sample_sex,sample_bam
C-48665L-N001-d,C-48665L,Normal,F,/path/to/C-48665L-N001-d.bam
C-PCYP90-N001-d,C-PCYP90,Normal,M,/path/to/C-PCYP90-N001-d.bam
C-MH6AL9-N001-d,C-MH6AL9,Normal,F,/path/to/C-MH6AL9-N001-d.bam
```

Here is an example command line usage for three samples:

```shell
biometrics extract \
  -i inputs.csv \
  --vcf /path/to/vcf \
  --bed /path/to/bed/file \
  -db /path/to/store/extract/output \
  -f /path/to/reference.fasta
```

#### Command line input
You can also specify each of your samples via command line arguments. Here is an example:

```shell
biometrics extract \
  -sn C-48665L-N001-d C-PCYP90-N001-d C-MH6AL9-N001-d \
  -sb /path/to/C-48665L-N001-d.bam /path/to/C-PCYP90-N001-d.bam /path/to/C-MH6AL9-N001-d.bam \
  -st Normal Normal Normal \
  -ss F M F \
  -sg C-48665L C-PCYP90 C-MH6AL9 \
  --vcf /path/to/vcf \
  --bed /path/to/bed/file \
  -db /path/to/store/extract/output \
  -f /path/to/reference.fasta
```
