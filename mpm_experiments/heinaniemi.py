import pandas as pd
import os
from src.allium_prepro.subtype_thesaurus import SubtypeThesaurus
from src.allium_prepro.gex_preprocessor import GexPreprocessor

dataset_name = 'heinaniemi'
print(f"Processing {dataset_name}...")
path_to_raw_data = '/home/mariya/Data/for_allium/raw/heinaniemi'
meta_input_file = f'{path_to_raw_data}/heinaniemi_meta.csv'
counts_input_file = \
    f'{path_to_raw_data}/GSE228632_RNAseq_read_counts.txt'

# Outputs
output_dir = '/home/mariya/Data/for_allium/allium'
pheno_output_file = f'{output_dir}/heinaniemi.pheno.allium.csv'
batches_output_file = f'{output_dir}/heinaniemi.batches.allium.csv'
counts_output_file = f'{output_dir}/heinaniemi.counts.raw.csv'

# PROCESS METADATA ########
print("Processing metadata...")

# Load the data
data = pd.read_csv(meta_input_file, index_col=0)

# Drop batch column and rename the rest
data = data.drop(columns=['batch', 'Subtype'])
data = data.rename(columns={'Subtype_updated': 'subtype',
                            'batches': 'batch'})
# Rename index
data.index.name = 'id'

# Extract batches info into a separate file
batches = data[['batch']]
batches.to_csv(batches_output_file, sep=',')

print("Processing phenotype data...")
# Drop batches column
data = data.drop(columns=['batch'])

# Get subtype translation
st = SubtypeThesaurus()

# Replace Subtype column using dict
data['subtype'] = st.translate_subtype_column(data['subtype'])

# Dump to output file
data.to_csv(pheno_output_file, sep=';')

# PROCESS COUNTS ########
# Load the data
data = pd.read_csv(counts_input_file, index_col=0, sep='\t')
data.to_csv(counts_output_file)

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='heinaniemi',
                    input_file=counts_output_file,
                    output_dir=output_dir,
                    gene_format='ensembl',
                    sample_col_regex='^(ALL|GE).*',
                    batches_file=batches_output_file)
p.run()
