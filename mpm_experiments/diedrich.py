
import pandas as pd
import os
from src.allium_prepro.subtype_thesaurus import SubtypeThesaurus
from src.allium_prepro.gex_preprocessor import GexPreprocessor

path_to_raw_data = '/home/mariya/Data/raw/diedrich'
pheno_input_file = f'{path_to_raw_data}/SupplementaryInfo.xlsx'
counts_input_file = \
    f'{path_to_raw_data}/GSE161501_ALL_cell_RNAseq_read_counts.txt'

# Outputs
output_dir = '/home/mariya/Data/allium'
pheno_output_file = f'{output_dir}/diedrich.pheno.allium.csv'
counts_output_file = f'{output_dir}/diedrich.counts.raw.csv'

# PROCESS PHENO ########
print("Processing phenotype data...")
# Load the data
data = pd.read_excel(pheno_input_file,
                     sheet_name='Supplemental Table 1',
                     header=2,
                     index_col='Sample SJ ID')

# Keep only index and subtype cols
data = data.iloc[:, :1]

# Rename columns
data = data.rename(columns={'Subtype_name': 'subtype'})
data.index.name = 'id'

# In subtype column, strip out everything after the underscore
data['subtype'] = data['subtype'].str.split('_').str[0]

# Get subtype translation
st = SubtypeThesaurus()
subtypes_dict = st.thesaurus()

# Replace Subtype column using dict
data['subtype'] = data['subtype'].replace(subtypes_dict)

# Dump to output file
data.to_csv(pheno_output_file, sep=';')

# PROCESS COUNTS ########
# We just need to change the separator from tab to comma and remove index name
data = pd.read_csv(counts_input_file, index_col=0, sep='\t')
data.index.name = None
data.to_csv(counts_output_file)

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='diedrich',
                    input_file=counts_output_file,
                    output_dir=output_dir,
                    gene_format='symbol',
                    sample_col_regex='^SJ.*')
p.run()

# Clean up intermediate file
os.remove(counts_output_file)
