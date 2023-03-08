---
description: For clustering samples into groups using genotype comparison data.
---

# Cluster

Takes as input the results from running `biometrics genotype` and clusters the samples together using the discordance rate. Done by thresholding the discordance rate into 0 or 1, where 1 means the sample pair came from the same patient. The default discordance rate threshold is 0.05, but you can change it via the `--discordance-threshold` argument. The tool then uses the `networkx` package to get the groups of samples that are connected.

{% hint style="info" %}
When you run `biometrics genotype`, it automatically outputs two sets of clustering results: (1) the first set just clusters your input samples, and (2) the second set clusters your input samples and samples in the database.
{% endhint %}

{% hint style="warning" %}
Due to the limitations of the discordance rate metric, samples that have contamination can lead to many false positive matches using this clustering approach. Hence, you might want to consider removing contaminated samples before running this tool.
{% endhint %}

## How to run the tool

You simply need to provide the CSV file that is outputted from running `biometrics genotype`.

```
biometrics cluster \
  -i genotype_comparison.csv \
  -o genotype_clusters.csv
```

You can also specify multiple inputs. It will automatically drop duplicate comparisons.

```
biometrics cluster \
  -i genotype_comparison_1.csv \
  -i genotype_comparison_2.csv \
  -i genotype_comparison_3.csv \
  -o genotype_clusters.csv
```

## Output

Produces a CSV file that contains the clustering results. Each row corresponds to a different sample. The table below provides a description on each column.

| Column Name                   | Description                                                                                                         |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| sample\_name                  | The sample name.                                                                                                    |
| expected\_sample\_group       | The expected group for the sample based on user input.                                                              |
| cluster\_index                | The integer cluster index. All rows with the same cluster\_index are in the same cluster.                           |
| cluster\_size                 | The size of the cluster this sample is in.                                                                          |
| avg\_discordance              | The average discordance between this sample and all other samples in the cluster.                                   |
| count\_expected\_matches      | The count of expected matches when comparing the sample to all others in the cluster.                               |
| count\_unexpected\_matches    | The count of unexpected matches when comparing the sample to all others in the cluster.                             |
| count\_expected\_mismatches   | The count of expected mismatches when comparing the sample to all other samples (inside and outside its cluster).   |
| count\_unexpected\_mismatches | The count of unexpected mismatches when comparing the sample to all other samples (inside and outside its cluster). |
