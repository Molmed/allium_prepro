
import pandas as pd
import os
from src.allium_prepro.batch_umap import BatchUmap

DATASET_PREFIX = 'heinaniemi'

# Data
data_dir = '/home/mariya/Data/allium'
pheno_file = f'{data_dir}/{DATASET_PREFIX}.pheno.allium.csv'
raw_counts_file = f'{data_dir}/{DATASET_PREFIX}.counts.raw.csv'
processed_counts_file = f'{data_dir}/{DATASET_PREFIX}.counts.allium.csv'
batches_file = f'{data_dir}/{DATASET_PREFIX}.batches.allium.csv'

bu = BatchUmap(prefix=f"{DATASET_PREFIX}_before",
               counts_file=raw_counts_file,
               batches_file=batches_file,
               phenotype_file=pheno_file,
               output_dir=data_dir,
               do_transform=True)
bu.run()

bu = BatchUmap(prefix=f"{DATASET_PREFIX}_after",
               counts_file=processed_counts_file,
               batches_file=batches_file,
               phenotype_file=pheno_file,
               output_dir=data_dir)
bu.run()
