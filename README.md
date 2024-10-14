# ALLIUM PrePro v1.0.0 :garlic:

## About

ALLIUM PrePro is a library for preprocessing gene expression (GEX) and DNA methylation (DNAm) data to prepare it for prediction using the [ALLIUM](https://github.com/Molmed/allium), a multimodal classifier of molecular subtypes in pediatric acute lymphoblastic leukemia.

### Publication

Krali, O., Marincevic-Zuniga, Y., Arvidsson, G. et al. Multimodal classification of molecular subtypes in pediatric acute lymphoblastic leukemia. npj Precis. Onc. 7, 131 (2023). https://doi.org/10.1038/s41698-023-00479-5

## Modules

This repository contains:
- GEX data preprocessing helpers
- metadata generation helpers (use only if changing reference genome versions)

DNAm preprocessing helpers are still in development.

## Conda environment
[Conda](https://docs.conda.io) must be installed on your system.

You will need to activate the `allium-prepro` conda environment before running any subsequent commands.

Install: `conda env create -f environment.yml`

Activate: `conda activate allium-prepro`

Update (after changes to environment.yml): `conda env update --file environment.yml --prune`

## Preprocessing GEX data
To prepare gene expression for prediction using ALLIUM, you will need a CSV file with raw gene transcript counts. The leftmost column should be HGNC gene symbols or Ensembl identifiers.

Example:
|         | Sample_1 | Sample_2 | ... |
| --------| -------- | -------- | --- |
| ETV6    | 10       | 10       | ... |
| SARS1   | 20       | 10       | ... |
| DOC2B   | 5        | 10       | ... |

### Batch correction

If your data was generated in batches, your count matrix should be batch corrected before further processing with ALLIUM PrePro. Alternatively, you may choose to pre-process the batches as separate count files and submit them separately to ALLIUM.

### Pre-processing and normalization for ALLIUM

Modify the `example_client.py` file in the project root, and run it with `python example_client.py`.

### Next steps
You are now ready to feed your PREFIX.counts.allium.csv file into [ALLIUM](https://github.com/Molmed/allium).
