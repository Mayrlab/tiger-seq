## Overview
This repository contains computational analysis used in Horste et al., 2023.

A preprint version of the article is [available on BioRxiv](https://doi.org/10.1101/2022.11.04.515216).

## Organization
The folders in the repository have the following purposes:

- `analysis` - primary source code and rendered HTMLs of R Markdown
- `data` - *input* data needed to run the code
- `envs` - Conda environment YAML files for recreating the execution environment
- `img` - *output* images
- `tbl` - *output* tables

### Source Code
The primary source code is found in the `analysis` folder. 
Folders and files are numbered in the expected order of execution.
The outputs of the `01-preprocess` files are already included in the `data` folder,
and therefore are not needed to be rerun. Anyone wishing to do further analysis or
replicate the models should be able to immediately run the `03-models` files.

### Input Data
The `data` folder contains the input data. We provide clean versions. Raw data versions
(needed to run the `01-preprocess` code) can be provide upon reasonable request.

### Execution Environment
The R instance used execute the files is captured both in the rendered RMDs themselves
(see **Runtime Details** section in HTMLs) and provided as YAML files in the `

To recreate on arbitrary platforms (Linux or MacOS), we recommend using [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#)
and the minimal YAML (`*.min.yaml`):

```
micromamba create -n brms_r41 -f envs/brms_r41.min.yaml
micromamba activate brms_r41
```

A full capture is also provided (`*.full.yaml`). This is only expected to recreate 
on the **osx-64** platform and is primarly intended for *exact* replication and a 
statement of record.

### Outputs
The `img` and `tbl` folders have output files that the source code will generate.
The `analysis` folder also has HTML renders of the R Markdown files.