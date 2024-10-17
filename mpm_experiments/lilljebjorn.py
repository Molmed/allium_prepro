
import pandas as pd
import os
from src.allium_prepro.subtype_thesaurus import SubtypeThesaurus
from src.allium_prepro.gex_preprocessor import GexPreprocessor

# GEX data from Lilljebjorn et al. 2016
path_to_raw_data = '/home/mariya/Data/raw/lilljebjorn'
pheno_input_file = f'{path_to_raw_data}/lilljebjorn.pheno.csv'
counts_input_file = f'{path_to_raw_data}/BCP-ALL_expected_counts.csv'

# Outputs
output_dir = '/home/mariya/Data/allium'
pheno_output_file = f'{output_dir}/lilljebjorn.pheno.allium.csv'
counts_output_file = f'{output_dir}/lilljebjorn.counts.raw.csv'

# PROCESS PHENO ########
print("Processing phenotype data...")
# Load the data
data = pd.read_csv(pheno_input_file, index_col=0, sep=';')

# Rename Subtype to subtype
data = data.rename(columns={'Subtype': 'subtype'})

# Rename index to public_id
data.index.name = 'id'

# For each integer index value, convert it to a string of the format "Case_00n"
data.index = data.index.map(lambda x: f'Case_{x:03d}')

# Get subtype translation
st = SubtypeThesaurus()
subtypes_dict = st.thesaurus()
data['subtype'] = st.translate_subtype_column(data['subtype'])

# Dumop to output file
data.to_csv(pheno_output_file, sep=';')

# PROCESS COUNTS ########
# Load the data
data = pd.read_csv(counts_input_file, index_col=0)

# Update index to use the gene name, after the first underscore
data.index = data.index.str.split('_', n=1).str[1]

# Drop all rows whose index starts with "ENSGR", these are Y chr genes
# https://www.gencodegenes.org/pages/faq.html
data = data[~data.index.str.startswith("ENSGR")]

# Write to file
data.to_csv(counts_output_file)

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='lilljebjorn',
                    input_file=counts_output_file,
                    output_dir=output_dir,
                    gene_format='symbol',
                    sample_col_regex='^Case.*')
p.run()

# Clean up intermediate file
os.remove(counts_output_file)
