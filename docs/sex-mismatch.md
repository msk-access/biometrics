---
description: Determine if a sample's predicted and known sex mismatch.
---

# Sex mismatch

This tool uses read coverage data on the Y chromosome to predict the sex for a sample, and then the compares the prediction to the expected sex to see if there is a mismatch. The metric requires the extracted coverage information from running the `extract` tool with the `--bed` flag supplied.

## How to run the tool

You can run this tool with one or more samples. At a minimum, the two required inputs are the list of sample names and the database \(biometrics will automatically load all sample data from the database\). Below is an example command:

```text
biometrics sexmismatch \
  -sn C-48665L-N001-d C-PCYP90-N001-d C-MH6AL9-N001-d \
  -db /path/to/store/extract/output
```

You can also indicate your input samples via a CSV file, which has the same format as what you provided to the extraction tool, but you only need the `sample_name` column:

```text
biometrics sexmismatch \
  -i samples.csv \
  -db /path/to/store/extract/output
```

## Output

All analyses output a CSV file containing the metrics for each sample. It will be saved either to the current working directory or to a folder you specify via `--outdir`. The table below describes each column in the CSV output.

| Column Name | Description |
| :--- | :--- |
| sample\_name | Sample name. |
| expected\_sex | The sample's expected sex. |
| predicted\_sex | The sample's predicted sex. |
| sex\_mismatch | True if expected and predicted sex mismatch. |

