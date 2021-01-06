---
description: Calculate minor contamination
---

# Minor contamination
Minor contamination is a metric to see if a sample is contaminated with small amounts of DNA from another unrelated sample. The metric requires the extracted pileup information.


### How to run the tool
You need one or more samples to run this analysis. The two required inputs are the list of sample names and the database (biometrics will automatically load all sample data from the database). Below is an example command:

```shell
biometrics minor \
  -sn C-48665L-N001-d C-PCYP90-N001-d C-MH6AL9-N001-d \
  -db /path/to/store/extract/output
```

You can also indicate your input samples via a CSV file, which has the same format as what you provided to the extraction tool, but you only need the `sample_name` column:

```shell
biometrics minor \
  -i samples.csv \
  -db /path/to/store/extract/output
```

### Output

All analyses output a CSV file containing the metrics for each sample. An interactive bar graph can also optionally be produced by supplying the `--plot` flag. These outputs are saved either to the current working directory or to a folder you specify via `--outdir`.

##### CSV file
The CSV file contains metrics for each pair of samples compared (one on each line). Table 1 below provides a description on each column.

| Column Name | Description |
| :--- | :--- |
| sample | Sample name. |
| sample_group | Sample group (if available). |
| sample_sex | Sample sex (if available). |
| sample_type | Sample type (if available). |
| total_homozygous_sites | Total number of homozygous sites. |
| minor_contamination | Minor contamination metric. |


##### Interactive plot
Below is an example bar plot showing the per-sample minor contamination metrics. You can hover over each bar to get more information about the sample. You can also control the minor contamination threshold (the horizontal red line) via the `--minor-threshold` flag.

![](.gitbook/assets/minor_contamination.html)

### Algorithm details

Minor contamination is calculated as the average minor allele frequency for homozygous sites. A homozygous site is defined as one with < 10% minor allele frequency.
