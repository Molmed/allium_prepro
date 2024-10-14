from lib.gex_preprocessor import GexPreprocessor

p = GexPreprocessor(prefix='jude',
                    input_file='/home/mariya/Data/jude/jude.counts.raw.csv',
                    output_dir='/home/mariya/Data/jude',
                    gene_format='symbol',
                    sample_col_regex='^SJ.*ALL.*')
p.run()
