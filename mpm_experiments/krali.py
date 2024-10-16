import pandas as pd
import os
from src.allium_prepro.subtype_thesaurus import SubtypeThesaurus
from src.allium_prepro.gex_preprocessor import GexPreprocessor

path_to_raw_data = '/home/mariya/Data/raw/krali'
meta_input_file = f'{path_to_raw_data}/41698_2023_479_MOESM2_ESM.xlsx'
counts_input_file = \
    f'{path_to_raw_data}/GSE227832_RNAseq_read_counts.txt'

# Outputs
output_dir = '/home/mariya/Data/allium'
pheno_output_file = f'{output_dir}/krali.pheno.allium.csv'
batches_output_file = f'{output_dir}/krali.batches.allium.csv'
counts_output_file = f'{output_dir}/krali.counts.raw.csv'

# PROCESS METADATA ########
print("Processing metadata...")

# Load the data
data = pd.read_excel(meta_input_file,
                     sheet_name='Supplementary Data 2',
                     header=1,
                     index_col='public_id')

# This column has a newline character that makes it tricky
subtype_col_name = data.columns[1]

data = data[[subtype_col_name,
             'GEX dataset',
             'library'
             ]]

# Rename cols
data = data.rename(columns={subtype_col_name: 'subtype',
                            'GEX dataset': 'partition',
                            'library': 'batch'})
# Rename index
data.index.name = 'id'

# Extract only hold out data
data = data[data['partition'] == 'held-out']

# Extract batches info into a separate file
batches = data[['batch']]
batches.to_csv(batches_output_file, sep=',')

print("Processing phenotype data...")
# Just keep subtype column
data = data[['subtype']]

# Get subtype translation
st = SubtypeThesaurus()
subtypes_dict = st.thesaurus()

# Replace Subtype column using dict
data['subtype'] = data['subtype'].replace(subtypes_dict)

# Dump to output file
data.to_csv(pheno_output_file, sep=';')

# Save subtypes data
subtypes = data

# FILTER FOR ONLY HOLDOUT DATA AND PROCESS COUNTS ########
# Load the data
data = pd.read_csv(counts_input_file, index_col=0, sep='\t')

# Only keep the columns that correlate to an id in subtypes
data = data[data.columns.intersection(subtypes.index)]

data.to_csv(counts_output_file)

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='krali',
                    input_file=counts_output_file,
                    output_dir=output_dir,
                    gene_format='ensembl',
                    sample_col_regex='^ALL.*',
                    batches_file=batches_output_file)
p.run()

# Clean up intermediate file
os.remove(counts_output_file)
