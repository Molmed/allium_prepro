import pandas as pd
import os
from src.allium_prepro.gex_concatenator import GexConcatenator
from src.allium_prepro.subtype_thesaurus import SubtypeThesaurus
from src.allium_prepro.gex_preprocessor import GexPreprocessor

dataset_name = 'tran'
print(f"Processing {dataset_name}...")
data_path = '/home/mariya/Data/raw/tran'
raw_data_dir = f'{data_path}/GSE181157_RAW/'
pheno_input_file = f'{data_path}/tran.pheno.csv'

processed_data_path = '/home/mariya/Data/allium'
pheno_output_file = f'{processed_data_path}/tran.pheno.allium.csv'
counts_output_file = f'{processed_data_path}/tran.counts.raw.csv'


# CONCAT GEX ########
def sample_name_extractor(x):
    return x.split('_', 1)[1].split('.')[0]


gc = GexConcatenator('tran',
                     raw_data_dir,
                     processed_data_path,
                     sample_name_extractor)
gc.concatenate()

# PROCESS PHENO ########
print("Processing phenotype data...")

# Load the data
data = pd.read_csv(pheno_input_file, index_col=0, sep=';')

# Rename subtype column
data = data.rename(columns={'Final subtype': 'subtype'})

# Get subtype translation
st = SubtypeThesaurus()
data['subtype'] = st.translate_subtype_column(data['subtype'])

# Rename index to public_id
data.index.name = 'id'

# Dumop to output file
data.to_csv(pheno_output_file, sep=';')

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='tran',
                    input_file=counts_output_file,
                    output_dir=processed_data_path,
                    gene_format='ensembl',
                    sample_col_regex='^16-.*')

p.run()
