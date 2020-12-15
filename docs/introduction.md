---
description: Basics on the usage of biometrics
---

# Introduction to Biometrics

Biometrics is a Python package to compute various metrics for assessing sample contamination, sample swaps, and sample sex validation. The package is composed of five tools:

#### Extract
Running this step is **required** before running any of the other four tools. This step preprocesses your BAM file(s) and stores the result in a file, which can then be accessed not just for your initial analysis but all other subsequent analyses making use of your samples. This provides a significant speed boost to running the four downstream biometrics tools.

Click [here](extraction.md) to read more about this tool.

#### Genotype
Compares each each sample against each other to verify expected sample matches and identify any unexpected matches or mismatches. Relies on computing a discordance score between each pair of samples.

Click [here](genotype.md) to read more about this tool.

#### Minor contamination
Minor contamination check is done to see if a patient’s sample is contaminated with little DNA from another unrelated individual.

Click [here](minor-contamination.md) to read more about this tool.

#### Major contamination
Major contamination check is done to see if a patient’s sample is contaminated with DNA from an unrelated individual.

Click [here](major-contamination.md) to read more about this tool.

#### Sex mismatch

Click [here](sex-mismatch.md) to read more about this tool.
