---
description: Basics on the usage of biometrics
---

# Introduction to Biometrics

Biometrics is a Python package to compute various metrics for assessing sample contamination, sample swaps, and sample sex validation. The package is composed of five tools (see below). All the tools (except the sex mismatch one) depend on you providing a VCF file of SNPs to use for computing the metrics. The sex mismatch tool requires you to provide a BED file containing the Y chromosome regions of interest.

#### Extract
Running this step is **required** before running any of the other four tools. This step extracts the pileup and coverage information from your BAM file(s) and stores the result in a file. The file can then be accessed not just for your initial analysis but for all subsequent analyses that make use of the sample. This provides a significant speed boost to running the four downstream biometrics tools.

Click [here](extraction.md) to read more about this tool.

#### Genotype
Compares each each sample against each other to verify expected sample matches and identify any unexpected matches or mismatches. Relies on computing a discordance score between each pair of samples.

Click [here](genotype.md) to read more about this tool.

#### Minor contamination
Minor contamination check is done to see if a patient’s sample is contaminated with a little DNA from unrelated individuals.

Click [here](minor-contamination.md) to read more about this tool.

#### Major contamination
Major contamination check is done to see if a patient’s sample is contaminated with DNA from unrelated individuals.

Click [here](major-contamination.md) to read more about this tool.

#### Sex mismatch
Used to determine if the predicted sex mismatches the expected sex for a given sample.

Click [here](sex-mismatch.md) to read more about this tool.
