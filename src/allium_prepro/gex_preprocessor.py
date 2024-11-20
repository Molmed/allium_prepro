import pandas as pd
import os
import re
from gene_thesaurus import GeneThesaurus
import rpy2.robjects as robjects


class GexPreprocessor():
    def __init__(self,
                 prefix,
                 input_file,
                 output_dir,
                 gene_format,
                 sample_col_regex,
                 batches_file=None,
                 ref_data_dir=None,
                 ref_genome='Homo_sapiens.GRCh38.103',
                 tmp_dir='/tmp'):

        self._input_file = input_file
        self._filtered_file_path = \
            f'{output_dir}/{prefix}.tmp.counts.filtered.csv'
        self._normalized_file_path = \
            f'{output_dir}/{prefix}.tmp.counts.norm.csv'
        self._output_file_path = f'{output_dir}/{prefix}.counts.allium.csv'
        self._missing_genes_path = f'{output_dir}/{prefix}.missing_genes.csv'

        # Optional batch correction
        self._batch_corrected_file_path = None
        self._batches_file_path = batches_file
        if self._batches_file_path:
            self._batch_corrected_file_path = \
                f'{output_dir}/{prefix}.tmp.counts.batch_corrected.csv'

        # Throw exception if gene_format is not 'symbol' or 'ensembl'
        if gene_format not in ['symbol', 'ensembl']:
            raise ValueError(
                'gene_format must be either "symbol" or "ensembl"')
        self._gene_format = gene_format
        self._sample_col_regex = sample_col_regex

        # If not ref data dir, use the local one
        if not ref_data_dir:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            ref_data_dir = os.path.join(script_dir, '../../data/reference')
        self._annot_file_path = os.path.join(
            ref_data_dir, (f'{ref_genome}.allium.annotations.filtered.csv'))

        self._gt = GeneThesaurus(data_dir=tmp_dir)

    def run(self):
        if self._batches_file_path:
            self._batch_correction()
        self._preprocess_genes()
        self._normalize()
        self._allium_format()
        self._cleanup()

    def _batch_correction(self):
        print('Correcting batch effects...')

        robjects.r('''
        library(edgeR)
        library(sva)

        batch_correct <- function(gex_path, batches_path, output_path) {
            x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)
            batches <- read.csv(batches_path, row.names = 1, header= TRUE, check.names = FALSE)
            # Ensure 'batches' is ordered to match 'x'
            batches <- batches[colnames(x), , drop = FALSE]
            batch = batches$batch
            correcteddata <- ComBat_seq(as.matrix(x), batch=batch)
            write.csv(correcteddata, output_path)
        }
        ''')

        r_batch_correct = robjects.r['batch_correct']
        r_batch_correct(self._input_file,
                        self._batches_file_path,
                        self._batch_corrected_file_path)

    def _preprocess_genes(self):
        print('Preprocessing genes...')

        input_file = self._input_file
        if self._batch_corrected_file_path:
            input_file = self._batch_corrected_file_path

        # Load the data
        data = pd.read_csv(input_file, index_col=0)
        ref = pd.read_csv(self._annot_file_path)

        if self._gene_format == 'ensembl':
            # Get gene names from ensembl ids
            translated_genes = self._gt.translate_genes(data.index.values,
                                                        source='ensembl_id',
                                                        target='symbol')

            def _standardize_ensembl_name(ensembl_name):
                if ensembl_name in translated_genes:
                    return translated_genes[ensembl_name]
                return ensembl_name

            data['gene_name_std'] = data.index.map(_standardize_ensembl_name)
        else:
            # Code to use if index is in symbol format
            # Get all gene name values to update to the latest standard
            updated_genes = self._gt.update_gene_symbols(data.index.values)

            # Copy data.index to gene_name_std
            data['gene_name_std'] = data.index

            # Update gene_name_std from updated_genes
            data['gene_name_std'] = data['gene_name_std'].map(
                updated_genes).fillna(data['gene_name_std'])

        # Do the same for ref
        updated_genes = self._gt.update_gene_symbols(ref['name'].values)

        # Set new column with the latest gene name
        ref['gene_name_std'] = ref['name'].map(
            updated_genes).fillna(ref['name'])

        filtered_rows = {}
        case_columns = [col for col in data.columns if
                        re.match(self._sample_col_regex, col)]

        def _find_ref_key(row):
            if row.name in ref['id'].values:
                return ref[ref['id'] == row.name]['id'].values[0]
            if row['gene_name_std'] in ref['gene_name_std'].values:
                ref_rows = ref[ref['gene_name_std'] == row['gene_name_std']]

                if len(ref_rows) == 1:
                    return ref_rows['id'].values[0]

                # See if one of the rows also matches a ref id
                for i, r in ref_rows.iterrows():
                    if r['id'] in ref['id'].values:
                        return r['id']

                # Return the first row
                return ref_rows['id'].values[0]

        for i, row in data.iterrows():
            # Strip whitespace
            key = _find_ref_key(row)

            # Print ref row length
            if key is not None:
                try:
                    key = key.strip()
                except Exception:
                    print(f'Error: {key}')
                    exit()

                # Just keep case cols in row
                row = row[case_columns]

                # Is there already a row in filtered_rows?
                if key in filtered_rows:
                    filtered_rows[key] += row
                else:
                    filtered_rows[key] = row
                # Keep the previous index name for the row
                filtered_rows[key].name = key

        data = pd.DataFrame(filtered_rows.values())

        # Print records in ref.id that are not in data.id
        missing = ref[~ref['id'].isin(data.index)]

        # Create records for all missing genes in data, filled with 0s
        # The missing$id is the index value, and all the case columns are 0
        missing_data = pd.DataFrame(index=missing['id'],
                                    columns=case_columns,
                                    data=0)

        # Dump missing data to file
        missing_data.to_csv(self._missing_genes_path)

        # Remove index name
        missing_data.index.name = None

        # Append the missing data to the data
        data = pd.concat([data, missing_data])

        # Sort data by index
        data = data.sort_index()

        # Dump data to file
        data.to_csv(self._filtered_file_path)

    def _normalize(self):
        print('Normalizing data...')

        robjects.r('''
        library(edgeR)
        library(sva)

        # create a function `get_cpm`
        normalize <- function(gex_path, ref_path, output_path) {
            x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)
            annot <- read.csv(ref_path, row.names = 1, header= TRUE, check.names = FALSE)

            x_length_norm <- ( (x*10^3 )/annot$length)
            d <- DGEList(counts=x_length_norm)
            TMM <- calcNormFactors(d, method="TMM")
            CPM <- cpm(TMM, log = TRUE)
            write.csv(CPM, output_path)
        }
        ''')

        r_normalize = robjects.r['normalize']
        r_normalize(self._filtered_file_path,
                    self._annot_file_path,
                    self._normalized_file_path)

    def _allium_format(self):
        print('Formatting data for ALLIUM...')

        # Load the data
        data = pd.read_csv(self._normalized_file_path, index_col=0)

        # Format
        data = data.T
        data.index.name = 'id'
        data.columns.name = None

        # Dump to file
        data.to_csv(self._output_file_path)

    def _cleanup(self):
        os.remove(self._filtered_file_path)
        os.remove(self._normalized_file_path)
        if self._batch_corrected_file_path:
            os.remove(self._batch_corrected_file_path)
        print('Done!')
