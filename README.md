# ALLIUM PrePro v1.0.0 :garlic:

## About

ALLIUM PrePro is a library for preprocessing gene expression (GEX) and DNA methylation (DNAm) data to prepare it for prediction using the [ALLIUM](https://github.com/Molmed/allium), a multimodal classifier of molecular subtypes in pediatric acute lymphoblastic leukemia.

### Publication

Krali, O., Marincevic-Zuniga, Y., Arvidsson, G. et al. Multimodal classification of molecular subtypes in pediatric acute lymphoblastic leukemia. npj Precis. Onc. 7, 131 (2023). https://doi.org/10.1038/s41698-023-00479-5

## Modules

This repository contains:
- GEX data preprocessing helpers
- metadata generation helpers (use only if ALLIUM has been re-trained with a different gene annotation version)

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

Modify `examples/example_gex_prepro.py`.
Run it with: `python -m examples.example_gex_prepro`.

### Next steps
You are now ready to feed your PREFIX.counts.allium.csv file into [ALLIUM](https://github.com/Molmed/allium).

### Reference pre-processing
Look at `examples/example_ref_prepro.py`. **Note!** The ReferencePreprocessor class only needs to be used in the event that ALLIUM has been re-trained using a different gene annotation version.

## MPM Experiments
Preprocessing for experiments in the [MPM Research Group](https://www.uu.se/en/department/medical-sciences/research/research-groups/molecular-precision-medicine) are in the `mpm_experiments` directory and can be replicated by running `python -m mpm_experiments.EXPERIMENT_NAME`
