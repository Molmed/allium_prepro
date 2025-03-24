from src.allium_prepro.gex_concatenator import GexConcatenator
from src.allium_prepro.gex_preprocessor import GexPreprocessor

# GEX CONCATENATION #
# If your GEX is a series of files, you can concatenate them into a single file
def sample_name_extractor(x):
    # Return everything before the first dot
    return x.split('.')[0]

gc = GexConcatenator('MYDATASET',
                     '/path/to/gex_files',
                     '/path/to/output',
                     sample_name_extractor)
gc.concatenate()

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
# sample_col_regex: regex to match the sample column names
p = GexPreprocessor(prefix='MYDATASET',
                    input_file='/path/to/MYDATASET.counts.raw.csv',
                    output_dir='/path/to/output',
                    gene_format='symbol',
                    sample_col_regex='^SAMPLENAME_PREFIX*',
                    batches_file='/path/to/batches.csv')
p.run()
