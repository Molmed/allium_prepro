from lib.gex_concatenator import GexConcatenator
from lib.gex_preprocessor import GexPreprocessor
from lib.jude_phenotype_parser import JudePhenotypeParser

# GEX CONCATENATION #
data_path = '/home/mariya/Data/raw/jude'
raw_data_dir = f'{data_path}/feature_counts'
raw_phenotype_path = f'{data_path}/SAMPLE_INFO.txt'
processed_data_path = '/home/mariya/Data/allium'


def sample_name_extractor(x):
    # Return everything before the first dot
    return x.split('.')[0]


def filename_filter_func(x):
    return not x.startswith('SJAML')


gc = GexConcatenator('jude',
                     raw_data_dir,
                     processed_data_path,
                     sample_name_extractor,
                     filename_filter_func=filename_filter_func)
gc.concatenate()

# PHENOTYPE PARSER ###
jpp = JudePhenotypeParser('jude', raw_phenotype_path, processed_data_path)
jpp.parse()

# GENE NAME STANDARDIZATION, COUNT NORMALIZATION, ALLIUM FORMATTING #
p = GexPreprocessor(prefix='jude',
                    input_file=f'{processed_data_path}/jude.counts.raw.csv',
                    output_dir=processed_data_path,
                    gene_format='symbol',
                    sample_col_regex='^SJ.*ALL.*')
p.run()
