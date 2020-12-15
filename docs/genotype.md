---
description: For finding sample matches and mismatches.
---

## Genotype
Running this step is a **prerequisite** before running any of the other tools. This preprocesses your BAM file(s) and stores the result into a python pickle file (which contains JSON data). This allows for much faster analyses that make repeated use of your samples. The output of this extraction step is stored in a folder that you specify with `-db` argument.
