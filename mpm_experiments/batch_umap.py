
import pandas as pd
import os
from src.allium_prepro.batch_umap import BatchUmap

DATASET_PREFIX = 'heinaniemi'

# Data
data_dir = '/home/mariya/Data/allium'
pheno_file = f'{data_dir}/{DATASET_PREFIX}.pheno.allium.csv'
counts_file = f'{data_dir}/{DATASET_PREFIX}.counts.allium.csv'
batches_file = f'{data_dir}/{DATASET_PREFIX}.batches.allium.csv'

bu = BatchUmap(prefix=DATASET_PREFIX,
               counts_file=counts_file,
               batches_file=batches_file,
               phenotype_file=pheno_file,
               output_dir=data_dir)
bu.run()
